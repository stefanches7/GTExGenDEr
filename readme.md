# GTExGenDEr

See [Google doc](https://docs.google.com/document/d/13my9tC9hNQITXJx0JcaoQhQj3OSR7DcZJSrnXeScHn0/edit?usp=sharing) for the project description

## Structure

This project uses [wBuild](https://wbuild.readthedocs.io) with [Snakemake](https://snakemake.readthedocs.io/en/stable/) for automatic pipeline execution. 
Build running `snakemake` CLI (ensure all the files required by `Snakefile` and `Scripts` are in place).
`Scripts` directory contains Scripts for the project, `Output/` HTML reports and processed data. `Output/html/index.html` is the HTML report entry point. From there, navigate to different analysis through the tabs on the header.

