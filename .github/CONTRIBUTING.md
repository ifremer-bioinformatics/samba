# samba: Contributing Guidelines

Hi and welcome!
Many thanks for taking an interest in improving samba.
Contributions to the code are even more welcome ;)

> If you need help using or modifying samba then the best place to ask is on our [assistance support]().

## Contribution workflow

If you'd like to write some code for samba, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [samba issues](https://github.com/ifremer-bioinformatics/samba/issues) to avoid duplicating work
    * If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [samba repository](https://github.com/ifremer-bioinformatics/samba) to your GitHub account
3. Make the necessary changes / additions within your forked repository following [Pipeline conventions](#pipeline-contribution-conventions)
4. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Tests

When you create a pull request with changes, be sure the samba tests still work.
Typically, pull-requests are only fully reviewed when the tests are passing,  though of course we can help out before then.

## Getting help

For further information/help, please consult the [samba documentation](https://github.com/ifremer-bioinformatics/samba/blob/master/docs/usage.md) and don't hesitate to get in touch on our [assistance support]().

## Pipeline contribution conventions

To make the samba code and processing logic more understandable for new contributors and to ensure quality, we semi-standardise the way the code and other contributions are written.

### Adding a new step

If you wish to contribute a new step, please use the following coding standards:

1. Define the corresponding input channel into your new process from the expected previous process channel
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new flags/options to all necessary `config files`.
5. Add any new flags/options to the help message
6. Add any new script in `bin/`.
7. Do local tests that the new code works properly and as expected.

### Default values

Parameters should be initialised / defined with default values all necessary `config files` under the `params` scope.

### Default processes resource requirements

Sensible defaults for process resource requirements (CPUs / memory / time) for a process should be defined in `conf/resources.config`. These should generally be specified generic with `withLabel:` selectors so they can be shared across multiple processes/steps of the pipeline. 

The process resources can be passed on to the tool dynamically within the process with the `${task.cpu}` and `${task.memory}` variables in the `script:` block.

### Naming schemes

Please use the following naming schemes, to make it easy to understand what is going where.

* initial process channel: `ch_output_from_<process>`
* intermediate and terminal channels: `ch_<previousprocess>_for_<nextprocess>`

### Software version reporting

If you add a new tool to the pipeline, please ensure you add the information of the tool to the `get_software_version` process.

Add to the script block of the process, something like the following:

```bash
<YOUR_TOOL> --version &> v_<YOUR_TOOL>.txt 2>&1 || true
```

or

```bash
<YOUR_TOOL> --help | head -n 1 &> v_<YOUR_TOOL>.txt 2>&1 || true
```
