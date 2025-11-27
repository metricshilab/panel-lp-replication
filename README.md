# Replication: Panel Local Projection


This repository hosts a replication package for

* Ziwei Mei, Liugang Sheng, Zhentao Shi (2025), "Nickell Bias in Panel Local Projection," *Journal of International Economics*. [[arxiv:2302.13455](https://arxiv.org/abs/2302.13455)].

* `paper_summary.pdf` is a minimum non-technical summary of the paper.

The replication code is written in R.

## Repository Structure and Replication Overview

- `applications/`: Empirical applications. Contains the data files (`empirical_*.csv`), the main R routine (`LP_panel_all.r`). The Jupyter notebook `replication.ipynb` reproduces the results and figures of all the four empirical application in the paper (Running time: <=30 seconds).
- `simulations/`: Monte Carlo simulations for the main text (Figures 1-3). See `simulations/README.md` for details about the scripts and figures.
- `simulations_appendix/`: Simulation code for Appendix C (Figures C1-C12). The workflow is explained in `simulations_appendix/readme.md` for instructions.

In short:

- Empirical results -> `applications/`
- Main-text simulations -> `simulations/`
- Appendix simulations -> `simulations_appendix/`


## R Package

A companion [R package](https://github.com/zhentaoshi/panel-local-projection) provides functions to implement the estimation method.

## Docker Environment

### DockerHub

An R environment is provided at [DockerHub](https://hub.docker.com/repository/docker/ztshi/plp). It can run in `Github`'s `codespaces` directly (Install a Jupyter Kernel extension in the remote vscode if popped up.)

<!-- Steps:

1. Open a cloud server, for example `github/codespaces`
2. Pull the image `docker pull ztshi/plp:v0.3`
3. Run the image. The port will be `8888`
4. In the web browser of jupyter interface, if login information is needed, copy the `string` after `http://127.0.0.1:8888/?token=<long_token_string>`
5. Done -->

### Local Environment
To run the R notebook in a fully reproducible environment, build the Docker image from the repository root:

```bash
docker build -t plp .
```

Launch the container, exposing JupyterLab on port `8888` and mounting the repository so that changes persist:

```bash
docker run --rm -it -p 8888:8888 \
  -e JUPYTER_TOKEN=panel-lp \
  -v ${PWD}:/home/jovyan/work \
  plp
```

The container will print a URL containing the token (or the `JUPYTER_TOKEN` you supplied). Open that link in a browser to access JupyterLab, then open `applications/replication.ipynb`. All required R packages (ggplot2, reshape2, ggpubr) are pre-installed in the image.

## Contributors

* [Pan Ji](https://github.com/PanJi-0), [Shen Shu](https://github.com/shenshuuu), [Shi Zhentao](https://github.com/zhentaoshi)
* Please contact [Shen Shu](https://github.com/shenshuuu) if you have any inquiries.


## License

This work is licensed under the MIT License.

