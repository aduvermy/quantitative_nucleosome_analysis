## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

To run nextflow on you computer you need to have java (>= 1.8) installed.

```sh
java --version
```

To be able to easily test tools already implemented for nextflow on your computer (`src/nf_modules/` to see their list). You need to have docker installed.

```sh
docker run hello-world
```

### Installing

To install nextflow on you computer simply run the following command:

```sh
src/install_nextflow.sh
```

Then to initialize a given tools run the following command:

```sh
src/docker_modules/<tool_name>/<tool_version>/docker_init.sh
```

For example to initialize `file_handle` version `0.1.1`, run:

```sh
src/docker_modules/file_handle/0.1.1/docker_init.sh
```

To initialize all the tools:
```sh
find src/docker_modules/ -name "docker_init.sh" | awk '{system($0)}'
```

## Running the tests

To run tests we first need to get a training set
```sh
cd data
git clone https://gitbio.ens-lyon.fr/LBMC/Hub/tiny_dataset.git
cp tiny_dataset/fastq/tiny_R1.fastq tiny_dataset/fastq/tiny2_R1.fastq
cp tiny_dataset/fastq/tiny_R2.fastq tiny_dataset/fastq/tiny2_R2.fastq
cp tiny_dataset/fastq/tiny_S.fastq tiny_dataset/fastq/tiny2_S.fastq
cd ..
```

Then to run the tests for a given tools run the following command:

```sh
src/nf_modules/<tool_name>/<tool_version>/tests.sh
```

For example to run the tests on `Bowtie2` run:

```sh
src/nf_modules/bowtie2/tests.sh
```
