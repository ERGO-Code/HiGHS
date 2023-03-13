Documentation for [HiGHS](https://github.com/ERGO-Code/HiGHS)

## Editing the Docs

To edit the documentation, checkout a branch and edit the Markdown files in
[docs/src].

To build locally, call
```
julia --project=. make.jl
```

and the website is generated in the `build/` folder. To check it out, load `build/index.html` in your browser.

When you are happy with the changes, rename the `build/` folder to `docs/` and push to the highs-docs branch. Alternatively, if you have pending changes you wish to discuss with the team, please checkout a new branch and open a PR to branch highs-docs.

The GH Page is generated from the `docs/` folder of the highs-docs branch.
