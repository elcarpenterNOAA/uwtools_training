"""
A driver for the GSI component.
"""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

from iotaa import asset, task, tasks

from uwtools.api.config import get_nml_config, get_yaml_config
from uwtools.config.formats.nml import NMLConfig
from uwtools.drivers.driver import DriverCycleBased
from uwtools.drivers.support import set_driver_docstring
from uwtools.strings import STR
from uwtools.utils.tasks import file, filecopy, symlink


class GSI(DriverCycleBased):
    """
    A driver for GSI.
    """

    # Workflow tasks

    @tasks
    def background_files(self):
        """
        The deterministic background estimate of the state.
        """
        bkg_config = self.config["background_files"]
        yield self.taskname("background files")
        paths = {}
        for i in range(bkg_config["io_layout"] + 1):
            fn_tmpls = get_yaml_config(bkg_config["files"])
            fn_tmpls.dereferene(
                context={
                  **self.config,
                  **self.cycle,
                  **{"iii": f{i:03d}},
                  }
                )
            paths.update(fn_tmpls)
        if bkg_config["protocol"] == "link":
            yield [
                symlink(target=t, linkname=self.rundir / l) for l, t in paths.items()
                ]
        else:
            yield [
                filecopy(src=s, dst=self.rundir / d) for d, s in paths.items()
                ]

    @task
    def coupler_res(self):
        """
        The coupler.res file.
        """
        fn = "coupler.res"
        yield self.taskname(fn)
        path = self.rundir / fn
        yield asset(path, path.is_file)
        template_file = Path(self.config[fn]["template_file"])
        yield file(template_file)
        render(
            input_file=template_file,
            output_file=path,
            overrides={
                **self.config[fn].get("template_values", {}),
                },
        )


    @tasks
    def crtm_files(self):
        """
        Common CRTM files.
        """
        yield self.taskname("common CRTM files")
        paths = self.config["crtm_files"]["files"]
        if self.config["crtm_files"]["protocol"] == "link":
            yield [
                symlink(target=t, linkname=self.rundir / l) for l, t in paths.items()
                ]
        else:
            yield [
                filecopy(src=s, dst=self.rundir / d) for d, s in paths.items()
                ]

    @tasks
    def crtm_satinfo_files(self):
        """
        CRTM files specific to the satinfo config file.
        """
        yield self.taskname("satinfo CRTM files")
        source = Path(self.config["crtm_satinfo_files"]["source"])
        satinfo_sats = self._get_sats.ref
        paths = {
                f"crtm_coeffs/{sat}.SpcCoeff.bin": source / f"{sat}.SpcCoeff.bin",
                f"crtm_coeffs/{sat}.TauCoeff.bin": source / f"{sat}.TauCoeff.bin"
                for sat in satinfo_sats
                }
        if self.config["crtm_satinfo_files"]["protocol"] == "link":
            yield [
                symlink(target=t, linkname=self.rundir / l) for l, t in paths.items()
                ]
        else:
            yield [
                filecopy(src=s, dst=self.rundir / d) for d, s in paths.items()
                ]


    @tasks
    def files_copied(self):
        """
        Files copied for run.
        """
        yield self.taskname("files copied")
        yield [
            filecopy(src=Path(src), dst=self.rundir / dst)
            for dst, src in self.config.get("files_to_copy", {}).items()
        ]

    @tasks
    def files_linked(self):
        """
        Files linked for run.
        """
        yield self.taskname("files linked")
        yield [
            symlink(target=Path(target), linkname=self.rundir / linkname)
            for linkname, target in self.config.get("files_to_link", {}).items()
        ]

    @task
    def global_ensemble_background_files(self):
        """
        The listed ensemble background files most recently available.
        """
        fn = "filelist03"
        yield self.taskname("global ensemble file list")
        path = self.rundir / fn
        yield asset(path, path.is_file)
        prior_cycles = self.config["global_ensemble_background_files"]["prior_cycles"]
        # Find the 9 hour forecast from those available that is closest in time to the current
        # cycle.
        # Likely running in real time, so look at the most recent time first
        ens_files = sorted(glob.glob(self.config["global_ensemble_background_files"]["ens_path"]),
                reverse=True)
        # Find the most recent synoptic time
        i = cycle.hour - cycle.hour % 6
        while i < prior_cycles:
            prev_cycle = cycle - timedelta(hours=i)
            ens = get_yaml_config(self.config["global_ensemble_background_files"]["path"])
            ens.dereference(
                context={
                    **self.config,
                    **{"prev_cycle": prev_cycle},
                    }
                )
            i++6

        yield None



    @task
    def namelist_file(self):
        """
        The namelist file.
        """
        path = self._input_config_path
        yield self.taskname(path.name)
        yield asset(path, path.is_file)
        base_file = self.config[STR.namelist].get(STR.basefile)
        # add a dependency on anything that might require a namelist change.
        # are there enough ensemble files? etc.
        yield file(Path(base_file)) if base_file else None
        self.create_user_updated_config(
            config_class=NMLConfig,
            config_values=self.config[STR.namelist],
            path=path,
            schema=self.namelist_schema(),
        )
        #TODO: tack on the stupid section

    @tasks
    def provisioned_rundir(self):
        """
        Run directory provisioned with all required content.
        """
        yield self.taskname("provisioned run directory")
        task_list = [
            self.background_files(),
            self.coupler_res(),
            self.crtm_files(),
            self.crtm_satinfo_files(),
            self.files_copied(),
            self.files_linked(),
            self.namelist_file(),
            self.observations(),
            self.provider_list(),
            self.radstat(),
            self.reject_list(),
            self.runscript(),
            self.satbias(),
            self.use_list(),
        ]
        if self.config["use_regional_ensemble"]:
            task_list.append(self.regional_ensemble_background_files())
        if self.config["use_global_ensemble"]:
            task_list.append(self.global_ensemble_background_files())
        yield task_list

    @task
    def observations(self):
        """
        The observation files.
        """
        yield self.taskname("observation files")
        paths = self.config["observations"]["files"]
        if self.config["observations"]["protocol"] == "link":
            yield [
                symlink(target=t, linkname=self.rundir / l) for l, t in paths.items()
                ]
        else:
            yield [
                filecopy(src=s, dst=self.rundir / d) for d, s in paths.items()
                ]

    @external
    def provider_list(self):
        """
        An optional provider list file.
        """
        fn = "gsd_sfcobs_provider.txt"
        yield self.taskname(fn)
        path = self.rundir / fn
        for src in self.config["provider_list"]:
            if (src_path := Path(src)).is_file():
                # Use the first one that exists
                yield asset(src_path, src_path.is_file())
        # It's okay if it's missing
        yield asset(None, lambda: True)

    @task
    def unzip(src: str, dst: Path, check: bool = True, clean: bool = False):
        """
        Unzip a file in the local filesystem.

        :param src: Path to the zipped file.
        :param dst: Path to the destination file to create.
        :param check: Check the existence of the zipped file before trying to unzip.
        :param clean: Remove src file after unziping?
        """
        yield "Local %s unzip -> %s" % (src, dst)
        yield asset(Path(dst), Path(dst).is_file)
        yield file(src) if check else None
        with gzip.open(src, mode="rt") as zipfile:
            Path(dst).write_text(zipfile.read())
        if clean:
            Path(src).unlink()

    @tasks
    def diag_files(self):
        """
        The first-guess diag files.
        """
        yield self.taskname("diag files")
        files = self.config.get("files", {})
        if files:
            zipfiles = {}
            regfiles = {}
            for dst, src in files.items():
                src_path = Path(src)
                if src_path.suffix == ".gz":
                    zipfiles[dst] = src
                else:
                    regfiles[dst] = src
            yield [
                filecopy(src=s, dst=self.rundir / d) for d, s in regfiles.items(),
                unzip(src=s, dst=self.rundir / d) for d, s in zipfiles.items(),
                ]
        if self.config.get("tarfile"):
            yield [self._untar_diag_files()]
        yield None

    @tasks
    def _untar_unzip_diag_files(tarfile: str):
        """
        The untarred first-guess diag files.
        """
        yield "Untar %s" % (tarfile)
        with tarfile.open(tarfile) as tf:
            filenames = tf.getmembers()
            tf.extractall(path=self.rundir)
        zipfiles = {fn.replace("_ges", ""): self.rundir / fn for fn in filenames}
        yield [unzip(src=s, dst=self.rundir / d, clean=True) for d, s in zipfiles.items()]

    @task
    def regional_ensemble_background_files(self):
        """
        The linked ensemble background files most recently available.
        """
        yield self.taskname("regional ensemble files")
        paths = self.config["regional_ensemble_background_files"]["files"]
        if self.config["regional_ensemble_background_files"]["protocol"] == "link":
            yield [
                symlink(target=t, linkname=self.rundir / l) for l, t in paths.items()
                ]
        else:
            yield [
                filecopy(src=s, dst=self.rundir / d) for d, s in paths.items()
                ]

    @task
    def reject_list(self):
        """
        An optional reject list file.
        """
        fn = "current_bad_aircraft"
        yield self.taskname(fn)
        path = self.rundir / fn
        for src in self.config["reject_list"]:
            if (src_path := Path(src)).is_file():
                # Use the first one that exists
                yield asset(src_path, src_path.is_file())
        # It's okay if it's missing
        yield asset(None, lambda: True)

    @task
    def satbias(self):
        """
        The most relevant satbias file.
        """

    @external
    def use_list(self):
        """
        An optional use list file.
        """
        fn = "gsd_sfcobs_uselist.txt"
        yield self.taskname(fn)
        path = self.rundir / fn
        for src in self.config["use_list"]:
            if (src_path := Path(src)).is_file():
                # Use the first one that exists
                log.info("Using {src} as use list")
                yield asset(src_path, src_path.is_file())
                return
        # It's okay if it's missing
        log.warning("No use list is available")
        yield asset(None, lambda: True)

# Consider a "go anyway" config option for the driver.

    # Public helper methods

    @classmethod
    def driver_name(cls) -> str:
        """
        The name of this driver.
        """
        return STR.gsi

    @property
    def output(self) -> dict[str, Path] | dict[str, list[Path]]:
        """
        Returns a description of the file(s) created when this component runs.
        """

    # Private helper methods

    @task
    def _get_sats(self):
        """
        Parse the satinfo file to get the list of satellites.
        """
        yield self.taskname("get sats from satinfo")
        satinfo = self.rundir / "satinfo"
        sats = ()
        yield asset(sats, lambda s: len(s))
        yield [self.files_to_link(), self.files_to_copy()]
        contents = satinfo.read_text().strip().split("\n")
        sats = sorted(set([c.split()[0] for c in contents if c.split()[0][0] != "!"]))


    @property
    def _input_config_path(self) -> Path:
        """
        Path to the input config file.
        """
        return self.rundir / "gsiparm.nml"


    @property
    def _runcmd(self) -> str:
        """
        The full command-line component invocation.
        """
        execution = self.config.get(STR.execution, {})
        mpiargs = execution.get(STR.mpiargs, [])
        components = [
            execution.get(STR.mpicmd),
            *[str(x) for x in mpiargs],
            "%s < %s" % (execution[STR.executable], self._input_config_path),
        ]
        return " ".join(filter(None, components))

set_driver_docstring(GSI)
