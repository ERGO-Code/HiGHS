# How to contribute to HiGHS

Welcome! This document explains some of the ways you can contribute to HiGHS.

## Code of Conduct

This project and everyone participating in it is governed by the [HiGHS Code of Conduct](https://github.com/ERGO-Code/HiGHS/blob/master/CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code.

## Contribute code to HiGHS

HiGHS is open source for distribution rather than contribution. This applies particularly to the core C++ code of the solvers. However, there is definitely scope for external contribution to interfaces and documentation. If you want to contribute in this way, please open an issue before making a pull request, since pull requests to the HiGHS solvers will not normally be accepted.

Please note that any pull requests to HiGHS must be made to the `latest` branch, so that full testing can take place before updating `master`. This ensures that `master` is always good to be downloaded by users. As a consequence, any changes should be developed as a branch of `latest` to reduce the chance of merge conflicts.

Note that, under the terms of the MIT license, by contributing to HiGHS you assign away your rights to the content of your contribution.

## AI generated contributions

The number of AI generated issues has increased dramatically since the summer of 2026. Relatively few will lead to practical contributions to the codebase, and a high volume of them obstructs the handling of those that do. If you want to submit AI generated issues, please:

- **Limit the number of issues.** Open no more than a few at a time. Large batches of issues, for example from an automated scan of the codebase, will not be accepted.
- **Prioritise.** Submit those that matter most to your use of HiGHS, and combine closely related points into a single issue.
- **Email us if you have more than a few issues.** Send a short list, in order of priority, to hello@highs.dev, and include your GitHub username so that we can recognise the source of the issues and invite their submission.
- **Check before submitting.** Search the open and closed issues for duplicates, confirm that the problem occurs with the latest version of HiGHS, and state in the issue that AI was used.

AI generated issues that do not follow this process may be closed or deleted without discussion.

Whilst AI may be valuable for generating and communicating suggestions for code and build system developments, AI-generated pull requests (other than the most trivial) will not be accepted. Reviewing such a pull request in detail takes longer than making the change directly, so if the suggestions are good, they will be implemented by HiGHS developers.

## Contact the HiGHS team

Contact [Julian](https://github.com/jajhall) (General issues and solvers), [Ivet](https://github.com/galabovaa) (Build and interfaces).

## Improve the documentation

The top level [documentation](https://ergo-code.github.io/HiGHS/) is created using [Documenter](https://documenter.juliadocs.org/stable/), with the files held on the [HiGHS repository](https://github.com/ERGO-Code/HiGHS/tree/master/docs). If your change is small (like fixing typos, or one or two sentence corrections), the easiest way to do this is to fork the repository, create a branch from `latest`, and then create a pull request. (See *Contribute code to HiGHS* for more on this.) If your change is larger, or touches multiple files, please raise an issue describing what you want to do.

## Raise an issue

You can raise an [issue](https://github.com/ERGO-Code/HiGHS/issues) with HiGHS in the normal way via [GitHub](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-an-issue). If AI was used to find or write up the issue, first read [AI generated contributions](#ai-generated-contributions).
