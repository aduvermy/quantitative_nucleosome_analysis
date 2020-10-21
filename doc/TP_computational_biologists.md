---
title: "TP for computational biologists"
author: Laurent Modolo [laurent.modolo@ens-lyon.fr](mailto:laurent.modolo@ens-lyon.fr)
date: 20 Jun 2018
output:
pdf_document:
toc: true
toc_depth: 3
    number_sections: true
highlight: tango
    latex_engine: xelatex
---

The goal of this practical is to learn how to *wrap* tools in [Docker](https://www.docker.com/what-docker) or [Environment Module](http://www.ens-lyon.fr/PSMN/doku.php?id=documentation:tools:modules) to make them available to nextflow on a personal computer or at the [PSMN](http://www.ens-lyon.fr/PSMN/doku.php).

Here we assume that you followed the [TP for experimental biologists](./TP_experimental_biologists.md), and that you know the basics of [Docker containers](https://www.docker.com/what-container) and [Environment Module](http://www.ens-lyon.fr/PSMN/doku.php?id=documentation:tools:modules). We are also going to assume that you know how to build and use a nextflow pipeline from the template [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow).

For the practical you can either work with the WebIDE of Gitlab, or locally as described in the [git: basis formation](https://gitlab.biologie.ens-lyon.fr/formations/git_basis).

# Docker

To run a tool within a [Docker container](https://www.docker.com/what-container) you need to write a `Dockerfile`.

[`Dockerfile`](./src/docker_modules/kallisto/0.44.0/Dockerfile) are found in the [pipelines/nextflow](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) project under `src/docker_modules/`. Each [`Dockerfile`](./src/docker_modules/kallisto/0.44.0/Dockerfile) is paired with a [`docker_init.sh`](./src/docker_modules/kallisto/0.44.0/docker_init.sh) file like following the example for `Kallisto` version `0.43.1`:

```sh
$ ls -l src/docker_modules/kallisto/0.43.1/
total 16K
drwxr-xr-x 2 laurent users 4.0K Jun 5 19:06 ./
drwxr-xr-x 3 laurent users 4.0K Jun 6 09:49 ../
-rw-r--r-- 1 laurent users  587 Jun  5 19:06 Dockerfile
-rwxr-xr-x 1 laurent users 79 Jun 5 19:06 docker_init.sh*
```

## [`docker_init.sh`](./src/docker_modules/kallisto/0.44.0/docker_init.sh)
The [`docker_init.sh`](./src/docker_modules/kallisto/0.44.0/docker_init.sh) is a simple sh script with executable rights (`chmod +x`). By executing this script, the user creates a [Docker container](https://www.docker.com/what-container) with the tool installed a specific version. You can check the [`docker_init.sh`](./src/docker_modules/kallisto/0.44.0/docker_init.sh) file of any implemented tools as a template.

Remember that the name of the [container](https://www.docker.com/what-container) must be in lower case and in the format `<tool_name>:<version>`.
For tools without a version number you can use a commit hash instead.

## [`Dockerfile`](./src/docker_modules/kallisto/0.44.0/Dockerfile)

The recipe to wrap your tool in a [Docker container](https://www.docker.com/what-container) is written in a [`Dockerfile`](./src/docker_modules/kallisto/0.44.0/Dockerfile) file.

For `Kallisto` version `0.44.0` the header of the `Dockerfile` is :

```Docker
FROM ubuntu:18.04
MAINTAINER Laurent Modolo

ENV KALLISTO_VERSION=0.44.0
```

The `FROM` instruction means that the [container](https://www.docker.com/what-container) is initialized from a bare installation of Ubuntu 18.04. You can check the versions of Ubuntu available [here](https://hub.docker.com/_/ubuntu/) or others operating systems like [debian](https://hub.docker.com/_/debian/) or [worst](https://hub.docker.com/r/microsoft/windowsservercore/).

Then we declare the *maintainer* of the container. Before declaring an environment variable for the container named `KALLISTO_VERSION`, which contains the version of the tool wrapped. This this bash variable will be declared for the user root within the [container](https://www.docker.com/what-container).

You should always declare a variable `TOOLSNAME_VERSION` that contains the version number of commit number of the tools you wrap. In simple cases you just have to modify this line to create a new `Dockerfile` for another version of the tool.

The following lines of the [`Dockerfile`](./src/docker_modules/kallisto/0.44.0/Dockerfile) are a succession of `bash` commands executed as the **root** user within the container.
Each `RUN` block is run sequentially by `Docker`. If there is an error or modifications in a `RUN` block, only this block and the following `RUN` will be executed.

You can learn more about the building of Docker containers [here](https://docs.docker.com/engine/reference/builder/#usage).

When you build your [`Dockerfile`](./src/docker_modules/kallisto/0.44.0/Dockerfile), instead of launching many times the [`docker_init.sh`](./src/docker_modules/kallisto/0.44.0/docker_init.sh) script to tests your [container](https://www.docker.com/what-container), you can connect to a base container in interactive mode to launch tests your commands.

```sh
docker run -it ubuntu:18.04 bash
KALLISTO_VERSION=0.44.0
```

# SGE / [PSMN](http://www.ens-lyon.fr/PSMN/doku.php)

To run easily tools on the [PSMN](http://www.ens-lyon.fr/PSMN/doku.php), you need to build your own [Environment Module](http://www.ens-lyon.fr/PSMN/doku.php?id=documentation:tools:modules).

You can read the Contributing guide for the [PMSN/modules](https://gitlab.biologie.ens-lyon.fr/PSMN/modules) project [here](https://gitlab.biologie.ens-lyon.fr/PSMN/modules/blob/master/CONTRIBUTING.md)

# Nextflow

The last step to wrap your tool is to make it available in nextflow. For this you need to create at least 4 files, like the following for Kallisto version `0.44.0`:

```sh
ls -lR src/nf_modules/kallisto
src/nf_modules/kallisto/:
total 12
-rw-r--r-- 1 laurent users 551 Jun 18 17:14 index.nf
-rw-r--r-- 1 laurent users 901 Jun 18 17:14 mapping_paired.nf
-rw-r--r-- 1 laurent users 1037 Jun 18 17:14 mapping_single.nf
-rwxr-xr-x 1 laurent users 627 Jun 18 17:14 tests.sh*
```

The [`.config` files](./src/nf_modules/kallisto/) file contains instructions for two profiles : `psmn` and `docker`.
The [`.nf` files](./src/nf_modules/kallisto/) file contains nextflow processes to use `Kallisto`.

The [`tests/tests.sh`](./src/nf_modules/kallisto/tests/tests.sh) script (with executable rights), contains a series of nextflow calls on the other `.nf` files of the folder. Those tests correspond to execution of the `*.nf` files present in the [`kallisto folder`](./src/nf_modules/kallisto/) on the [LBMC/tiny_dataset](https://gitlab.biologie.ens-lyon.fr/LBMC/tiny_dataset) dataset with the `docker` profile. You can read the *Running the tests* section of the [README.md](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/blob/master/README.md).

## [`kallisto.config`](./src/nf_modules/kallisto/)

The `.config` file defines the configuration to apply to your process conditionally to the value of the `-profile` option. You must define configuration for at least the `psmn` and `docker` profile.

```Groovy
profiles {
  docker {
    docker.temp = 'auto'
    docker.enabled = true
    process {
    }
  }
  psmn {
    process{
    }
  }
```

### `docker` profile

The `docker` profile starts by enabling docker for the whole pipeline. After that you only have to define the container name for each process:
For example, for `Kallisto` with the version `0.44.0`, we have:

```Groovy
process {
  withName: index_fasta {
    container = "kallisto:0.44.0"
  }
  withName: mapping_fastq {
    container = "kallisto:0.44.0"
  }
}
```

### `psmn` profile

The `psmn` profile defines for each process all the informations necessary to launch your process on a given queue with SGE at the [PSMN](http://www.ens-lyon.fr/PSMN/doku.php).
For example, for `Kallisto`, we have:

```Groovy
process{
  withName: index_fasta {
    beforeScript = "source /usr/share/lmod/lmod/init/bash; module use ~/privatemodules"
    module = "Kallisto/0.44.0"
    executor = "sge"
    cpus = 16
    memory = "30GB"
    time = "24h"
    queue = 'E5-2670deb128A,E5-2670deb128B,E5-2670deb128C,E5-2670deb128D,E5-2670deb128E,E5-2670deb128F'
    penv = 'openmp16'
  }
  withName: mapping_fastq {
    beforeScript = "source /usr/share/lmod/lmod/init/bash; module use ~/privatemodules"
    module = "Kallisto/0.44.0"
    executor = "sge"
    cpus = 16
    memory = "30GB"
    time = "24h"
    queue = 'E5-2670deb128A,E5-2670deb128B,E5-2670deb128C,E5-2670deb128D,E5-2670deb128E,E5-2670deb128F'
    penv = 'openmp16'
  }
}
```

The `beforeScript` variable is executed before the main script for the corresponding process.

## [`kallisto.nf`](./src/nf_modules/kallisto/kallisto.nf)

The [`kallisto.nf`](./src/nf_modules/kallisto/kallisto.nf) file contains examples of nextflow process that execute Kallisto.

- Each example must be usable as it is to be incorporated in a nextflow pipeline.
- You need to define, default value for the parameters passed to the process. 
- Input and output must be clearly defined.
- Your process should be usable as a starting process or a process retrieving the output of another process.

For more informations on processes and channels you can check the [nextflow documentation](https://www.nextflow.io/docs/latest/index.html).

## Making your wrapper available to the LBMC

To make your module available to the LBMC you must have a `tests.sh` script and one or many `docker_init.sh` scripts working without errors.
All the processes in your `.nf` must be covered by the tests.

After pushing your modifications on your forked repository, you can make a Merge Request to the [PSMN/modules](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow) **dev** branch. Where it will be tested and integrated to the **master** branch.

You can read more on this process [here](https://guides.github.com/introduction/flow/)


### `docker_modules`

The `src/docker_modules` contains the code to wrap tools in [Docker](https://www.docker.com/what-docker). [Docker](https://www.docker.com/what-docker) is a framework that allows you to execute software within [containers](https://www.docker.com/what-container). The `docker_modules` contains directory corresponding to tools and subdirectories corresponding to their version.

```sh
ls -l src/docker_modules/
rwxr-xr-x  3 laurent _lpoperator   96 May 25 15:42 bedtools/
drwxr-xr-x  4 laurent _lpoperator  128 Jun 5 16:14 bowtie2/
drwxr-xr-x  3 laurent _lpoperator   96 May 25 15:42 fastqc/
drwxr-xr-x  4 laurent _lpoperator  128 Jun 5 16:14 htseq/
```

To each `tools/version` corresponds two files:

```sh
ls -l src/docker_modules/bowtie2/2.3.4.1/
-rw-r--r-- 1 laurent _lpoperator  283 Jun  5 15:07 Dockerfile
-rwxr-xr-x  1 laurent _lpoperator   79 Jun 5 16:18 docker_init.sh*
```

The `Dockerfile` is the [Docker](https://www.docker.com/what-docker) recipe to create a [container](https://www.docker.com/what-container) containing `Bowtie2` in its `2.3.4.1` version. And the `docker_init.sh` file is a small script to create the [container](https://www.docker.com/what-container) from this recipe.

By running this script you will be able to easily install tools in different versions on your personal computer and use it in your pipeline. Some of the advantages are:

- Whatever the computer, the installation and the results will be the same
- You can keep [container](https://www.docker.com/what-container) for old version of tools and run it on new systems (science = reproducibility)
- You don’t have to bother with tedious installation procedures, somebody else already did the job and wrote a `Dockerfile`.
- You can easily keep [containers](https://www.docker.com/what-container) for different version of the same tools.

### `psmn_modules`

The `src/psmn_modules` folder is not really there. It’s a submodule of the project [PSMN/modules](https://gitlab.biologie.ens-lyon.fr/PSMN/modules). To populate it locally you can use the following command:

```sh
git submodule init
```

Like the `src/docker_modules` the [PSMN/modules](https://gitlab.biologie.ens-lyon.fr/PSMN/modules) project describe recipes to install tools and use them. The main difference is that you cannot use [Docker](https://www.docker.com/what-docker) on the PSMN. Instead you have to use another framework [Environment Module](http://www.ens-lyon.fr/PSMN/doku.php?id=documentation:tools:modules) which allows you to load modules for specific tools and version.
The [README.md](https://gitlab.biologie.ens-lyon.fr/PSMN/modules/blob/master/README.md) file of the [PSMN/modules](https://gitlab.biologie.ens-lyon.fr/PSMN/modules) repository contains all the instruction to be able to load the modules maintained by the LBMC and present in the [PSMN/modules](https://gitlab.biologie.ens-lyon.fr/PSMN/modules) repository.

## Create your Docker containers

For this practical, we are going to need the following tools:

- For Illumina adaptor removal: cutadapt
- For reads trimming by quality: UrQt
- For mapping and quantifying reads: BEDtools and Kallisto

To initialize these tools, follow the **Installing** section of the [README.md](https://gitlab.biologie.ens-lyon.fr/pipelines/nextflow/blob/master/README.md) file.

**If you are using a CBP computer don’t forget to clean up your docker containers at the end of the practical with the following commands:**

```sh
docker rm $(docker stop $(docker ps -aq))
docker rmi $(docker images -qf "dangling=true")
```


