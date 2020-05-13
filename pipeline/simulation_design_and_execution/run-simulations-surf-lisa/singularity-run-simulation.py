#!/usr/bin/env python

from subprocess import Popen
from multiprocessing import Pool
import argparse
import os

if __name__ == "__main__":
    # Get arguments
    parser = argparse.ArgumentParser(description='Run one or more singularity containers for my simulation study <https://github.com/JasperHG90/sleepsimR-run>')
    # Add the arguments
    parser.add_argument('host',
                        type=str,
                        help='Host on which the sleepsimR-API <https://github.com/JasperHG90/sleepsimR-api> is running')
    parser.add_argument('username',
                        type=str,
                        help="Username used to authenticate with the sleepsimR-API.")
    parser.add_argument('password',
                        type=str,
                        help="Password used to authenticate with the sleepsimR-API.")
    parser.add_argument('cmd',
                        type=str,
                        help="Command to execute.")
    parser.add_argument('--iterations',
                        type=int,
                        default=1,
                        help='Number of concurrent simulation iterations to run.')
    # Get arguments
    args = parser.parse_args()
    host = args.host
    usr = args.username
    pwd = args.password
    cmd = args.cmd
    iterations = args.iterations
    # Set environment variables
    os.environ["SINGULARITYENV_SLEEPSIMR_MASTER_HOST"] = host
    os.environ["SINGULARITYENV_SLEEPSIMR_API_PASSWORD"] = pwd
    os.environ["SINGULARITYENV_SLEEPSIMR_API_USERNAME"] = usr
    # Make several containers
    # Add sleep command for each 
    # This helps spread out the requests to the API
    slp_multiple = 1
    containers = []
    for idx in range(0, iterations):
        if idx == 0:
            containers.append(cmd)
        else:
            containers.append("sleep {} && {}".format(idx * slp_multiple, cmd))
    # Run processes
    processes = [Popen(cmd, shell=True) for cmd in containers]
    # Wait until complete
    for p in processes: 
        p.wait()
