# nextflow pipeline

This repository is a template and a library repository to help you build nextflow pipeline.
You can fork this repository to build your own pipeline.
To get the last commits from this repository into your fork use the following commands:

```sh
git remote add upstream gitlab_lbmc:pipelines/nextflow.git
git pull upstream master
```
**If you created your `.config` file before version `0.4.0` you need to run the script `src/.update_config.sh` to use the latest docker, singularity and conda configuration (don't forget to check your config files afterward for typos).**

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

[you can follow them here.](doc/getting_started.md)

## Available tools

[The list of available tools.](doc/available_tools.md)

## Projects using nextflow

[A list of projects using nextflow at the LBMC.](doc/nf_projects.md)

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://gitbio.ens-lyon.fr/pipelines/nextflow/tags).

## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitbio.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details
