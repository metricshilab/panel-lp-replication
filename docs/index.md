---
layout: default
title: "Overview"
nav_order: 1
permalink: /
---

# Panel Local Projection Replication
{: .fs-9 }

Ensuring transparent and reproducible results for *"Nickell Bias in Panel Local Projection: Financial Crises Are Worse Than You Think."*
{: .lead }

> Ziwei Mei, Liugang Sheng, Zhentao Shi (2025), *Journal of International Economics*
> [[arxiv:2302.13455](https://arxiv.org/abs/2302.13455)]

## What You'll Find Here

- **[Paper Summary](https://github.com/metricshilab/panel-lp-replication/blob/main/paper_summary.pdf)** - Non-technical minimum summary of the paper, with recommended procedure for applications

- **Empirical Examples** - R scripts in Jupyter notebook that replicate every empirical result in seconds.
- **Simulations** - Monte Carlo experiments behind the paper's main and appendix figures.
- **Docker workflow** - A hardened reproduction environment with every dependency baked in.

[Explore the Repository Contents]({{ "/contents/" | relative_url }}){: .btn .btn-primary .mr-2 }
[Set Up the Docker Environment]({{ "/docker/" | relative_url }}){: .btn }

## Why Local Projections Matter

The package implements the panel local projection estimator proposed in the paper, providing tools to study dynamic responses to shocks while addressing Nickell bias. With minimal setup you can rerun the authors' analysis.


# External Packages

Enhance or extend the analysis using the official packages:

| Platform | Repository | Notes |
| --- | --- | --- |
| R | [`panel-local-projection`](https://github.com/zhentaoshi/panel-local-projection) | Functions to estimate the panel LP model used in the paper. |
| Stata (in progress) | [`panel-local-projection-stata`](https://github.com/shenshuuu/panel-local-projection-stata) | Ongoing port of the estimator to Stata. |

These packages are optional for replication but useful when applying the method to new data sets.


## Maintainers

- [Pan Ji](https://github.com/PanJi-0)
- [Shen Shu](https://github.com/shenshuuu) - main contact
- [Shi Zhentao](https://github.com/zhentaoshi)
