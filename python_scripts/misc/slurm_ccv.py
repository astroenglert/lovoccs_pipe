
from lsst.ctrl.bps.parsl.sites.slurm import Slurm, Kwargs

from parsl.launchers import SrunLauncher

from parsl.executors import HighThroughputExecutor
from parsl.executors.base import ParslExecutor
from parsl.providers import SlurmProvider

from lsst.ctrl.bps.parsl.configuration import get_bps_config_value, get_workflow_name

# define a new class which inherits from Slurm with updated configs for use w. CCV
class SlurmCCV(Slurm):
    """
    
    Configuration for Slurm at Brown University's Center for Computation and Visualization
    
    This is, in-fact, identical to the base Slurm class, but with make_executor() overwritten 
    to disable the request for exclusive nodes and to use a different executor. Unfortunately,
    this does require manually overriding make_executor since the variables set in the executor
    and the provider are parameters, rather than attributes, and aren't accessible outward.
    
    """ 
    
    def make_executor(
        self,
        label: str,
        *,
        nodes: int | None = None,
        cores_per_node: int | None = None,
        walltime: str | None = None,
        mem_per_node: int | None = None,
        mem_per_worker: float | None = None,
        qos: str | None = None,
        constraint: str | None = None,
        singleton: bool = False,
        scheduler_options: str | None = None,
        provider_options: Kwargs | None = None,
        executor_options: Kwargs | None = None,
    ) -> ParslExecutor:
        """Return an executor for running on a Slurm cluster

        Parameters
        ----------
        label : `str`
            Label for executor.
        nodes : `int`, optional
            Default number of nodes for each Slurm job.
        cores_per_node : `int`, optional
            Default number of cores per node for each Slurm job.
        walltime : `str`, optional
            Default time limit for each Slurm job.
        mem_per_node : `float`, optional
            Memory per node (GB) to request for each Slurm job.
        mem_per_worker : `float`, optional
            Minimum memory per worker (GB), limited by the executor.
        qos : `str`, optional
            Quality of service for each Slurm job.
        constraint : `str`, optional
            Node feature(s) to require for each Slurm job.
        singleton : `bool`, optional
            Allow only a single Slurm job to run at a time?
        scheduler_options : `str`, optional
            ``#SBATCH`` directives to prepend to the Slurm submission script.
        provider_options : `dict`, optional
            Additional arguments for `SlurmProvider` constructor.
        executor_options : `dict`, optional
            Additional arguments for `HighThroughputExecutor` constructor.

        Returns
        -------
        executor : `HighThroughputExecutor`
            Executor for Slurm jobs.
        """
        nodes = get_bps_config_value(self.site, "nodes", int, nodes, required=True)
        cores_per_node = get_bps_config_value(self.site, "cores_per_node", int, cores_per_node)
        walltime = get_bps_config_value(self.site, "walltime", str, walltime, required=True)
        mem_per_node = get_bps_config_value(self.site, "mem_per_node", int, mem_per_node)
        qos = get_bps_config_value(self.site, "qos", str, qos)
        singleton = get_bps_config_value(self.site, "singleton", bool, singleton)
        scheduler_options = get_bps_config_value(self.site, "scheduler_options", str, scheduler_options)

        job_name = get_workflow_name(self.config)
        if scheduler_options is None:
            scheduler_options = ""
        else:
            scheduler_options += "\n"
        scheduler_options += f"#SBATCH --job-name={job_name}\n"
        if qos:
            scheduler_options += f"#SBATCH --qos={qos}\n"
        if constraint:
            scheduler_options += f"#SBATCH --constraint={constraint}\n"
        if singleton:
            # The following SBATCH directives allow only a single slurm job
            # (parsl block) with our job_name to run at once. This means we can
            # have one job running, and one already in the queue when the first
            # exceeds the walltime limit. More backups could be achieved with a
            # larger value of max_blocks. This only allows one job to be
            # actively running at once, so that needs to be sized appropriately
            # by the user.
            scheduler_options += "#SBATCH --dependency=singleton\n"
        return HighThroughputExecutor(
            label,
            provider=SlurmProvider(
                nodes_per_block=nodes,
                cores_per_node=cores_per_node,
                mem_per_node=mem_per_node,
                walltime=walltime,
                scheduler_options=scheduler_options,
                exclusive=False, # changed from parent class
                launcher=SrunLauncher(overrides="--export=ALL"), # changed from parent class
                **(provider_options or {}),
            ),
            mem_per_worker=mem_per_worker,
            address=self.get_address(),
            **(executor_options or {}),
        )
    
    
    
    

    
