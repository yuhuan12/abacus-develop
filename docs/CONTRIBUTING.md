# Contributing to ABACUS

First of all, thank you for taking time to make contributions to ABACUS!
This file provides the more technical guidelines on how to realize it.
For more non-technical aspects, please refer to the [ABACUS Contribution Guide](./community/contribution_guide.md)

## Table of Contents

- [Got a question?](#got-a-question)
- [Structure of the package](#structure-of-the-package)
- [Submitting an Issue](#submitting-an-issue)
- [Comment style for documentation](#comment-style-for-documentation)
- [Code formatting style](#code-formatting-style)
- [Generating code coverage report](#generating-code-coverage-report)
- [Adding a unit test](#adding-a-unit-test)
- [Debugging the codes](#debugging-the-codes)
- [Submitting a Pull Request](#submitting-a-pull-request)
- [Commit message guidelines](#commit-message-guidelines)

## Got a question?

Please referring to our GitHub [issue tracker](https://github.com/deepmodeling/abacus-develop/issues), and our developers are willing to help.
If you find a bug, you can help us by submitting an issue to our GitHub Repository. Even better, you can submit a Pull Request with a patch. You can request a new feature by submitting an issue to our GitHub Repository.
If you would like to implement a new feature, please submit an issue with a proposal for your work first, and that ensures your work collaborates with our development road map well. For a major feature, first open an issue and outline your proposal so that it can be discussed. This will also allow us to better coordinate our efforts, prevent duplication of work, and help you to craft the change so that it is successfully accepted into the project.

## Structure of the package

Please refer to [our instructions](./quick_start/easy_install.md) on how to installing ABACUS.
The source code of ABACUS is based on several modules. Under the ABACUS root directory, there are the following folders:

- `cmake`: relevant files for finding required packages when compiling the code with cmake;
- `docs`: documents and supplementary info about ABACUS;
- `examples`: some examples showing the usage of ABACUS;
- `source`: the source code in separated modules, under which a `test` folder for its unit tests;
- `tests`: End-to-end test cases;
- `tools`: the script for generating the numerical atomic orbitals.

## Submitting an Issue

Before you submit an issue, please search the issue tracker, and maybe your problem has been discussed and fixed. You can [submit new issues](https://github.com/deepmodeling/abacus-develop/issues/new/choose) by filling our issue forms.
To help us reproduce and confirm a bug, please provide a test case and building environment in your issue.

## Comment style for documentation

ABACUS uses Doxygen to generate docs directly from `.h` and `.cpp` code files.

For comments that need to be shown in documents, these formats should be used -- **Javadoc style** (as follow) is recommended, though Qt style is also ok. See it in [official manual](https://www.doxygen.nl/manual/docblocks.html).

A helpful VS Code extension -- [Doxygen Documentation Generator](https://marketplace.visualstudio.com/items?itemName=cschlosser.doxdocgen), can help you formating comments.

An practical example is class [LCAO_Deepks](https://github.com/deepmodeling/abacus-develop/blob/deepks/source/module_hamilt_lcao/module_deepks/LCAO_deepks.h), the effects can be seen on [readthedocs page](https://abacus-deepks.readthedocs.io/en/latest/DeePKS_API/classLCAO__Descriptor.html#exhale-class-classLCAO-Descriptor)

- Tips
  - Only comments in .h file will be visible in generated  by Doxygen + Sphinx;
  - Private class members will not be documented;
  - Use [Markdown features](https://www.doxygen.nl/manual/markdown.html), such as using a empty new line for a new paragraph.

- Detailed Comment Block

    ```cpp
    /**
    * ... text ...
    */
    ```

- Brief + Detailed Comment Block

    ```cpp
    /// Brief description which ends at this dot. Details follow
    /// here.

    /// Brief description.
    /** Detailed description. */
    ```

- Comments After the Item: Add a "<"

    ```cpp
    int var; /**<Detailed description after the member */
    int var; ///<Brief description after the member
    ```

- Parameters
    usage: `[in],[out],[in,out] description`
    *e.g.*

    ```cpp
    void foo(int v/**< [in] docs for input parameter v.*/);
    ```

    or use `@param` command.

- Formula
  - inline: `\f$myformula\f$`
  - separate line: `\f[myformula\f]`
  - environment: `\f{environment}{myformula}`
  - *e.g.*

    ```latex
    \f{eqnarray*}{
            g &=& \frac{Gm_2}{r^2} \\
            &=& \frac{(6.673 \times 10^{-11}\,\mbox{m}^3\,\mbox{kg}^{-1}\,
                \mbox{s}^{-2})(5.9736 \times 10^{24}\,\mbox{kg})}{(6371.01\,\mbox{km})^2} \\
            &=& 9.82066032\,\mbox{m/s}^2
    \f}
    ```

## Code formatting style

We use `clang-format` as our code formatter. The `.clang-format` file in root directory describes the rules to conform with.
For Visual Studio Code developers, the [official extension of C/C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) provided by Microsoft can help you format your codes following the rules. With this extension installed, format your code with `shift+command/alt+f`.
Configure your VS Code settings as `"C_Cpp.clang_format_style": "file"` (you can look up this option by pasting it into the search box of VS Code settings page), and all this stuff will take into effect. You may also set `"editor.formatOnSave": true` to avoid formatting files everytime manually.

## Adding a unit test

We use GoogleTest as our test framework. Write your test under the corresponding module folder at `abacus-develop/tests`, then append the test to `tests/CMakeLists.txt`. If there are currently no unit tests provided for the module, do as follows. `module_base` provides a simple demonstration.

- Add a folder named `test` under the module.
- Append the content below to `CMakeLists.txt` of the module:

    ```cmake
    IF (BUILD_TESTING)
      add_subdirectory(test)
    endif()
    ```

- Add a blank `CMakeLists.txt` under `module*/test`.

To add a unit test:

- Write your test under `GoogleTest` framework.
- Add your testing source code with suffix `*_test.cpp` in `test` directory.
- Append the content below to `CMakeLists.txt` of the module:

    ```cmake
    AddTest(
      TARGET <module_name>_<test_name> # this is the executable file name of the test
      SOURCES <test_name>.cpp

      # OPTIONAL: if this test requires external libraries, add them with "LIBS" statement.
      LIBS math_libs # `math_libs` includes all math libraries in ABACUS.
    )
    ```

- Build with `-D BUILD_TESTING=1` flag. You can find built testing programs under `build/source/<module_name>/test`.
- Follow the installing procedure of CMake. The tests will move to `build/test`.

## Debugging the codes

For the unexpected results when developing ABACUS, [GDB](https://www.sourceware.org/gdb/) will come in handy.

1. Compile ABACUS with debug mode.

    ```bash
    cmake -B build -DCMAKE_BUILD_TYPE=Debug
    ```

2. After building and installing the executable, enter the input directory, and launch the debug session with `gdb abacus`. For [debugging in Visual Studio Code](https://code.visualstudio.com/docs/cpp/cpp-debug), please set [cwd](https://code.visualstudio.com/docs/cpp/launch-json-reference#_cwd) to the input directory, and [program](https://code.visualstudio.com/docs/cpp/launch-json-reference#_program-required) to the path of ABACUS executable.

3. Set breakpoints, and run ABACUS by typing "run" in GDB command line interface. If the program hits the breakpoints or exception is throwed, GDB will stop at the erroneous code line. Type "where" to show the stack backtrace, and "print i" to get the value of variable i.

4. For debugging ABACUS in multiprocessing situation, `mpirun -n 1 gdb abacus : -n 3 abacus` will attach GDB to the master process, and launch 3 other MPI processes.

For segmentation faults, ABACUS can be built with [Address Sanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer) to locate the bugs. This feature requires a GCC or Clang compiler, and does not support Intel compiler.

```bash
cmake -B build -DENABLE_ASAN=1
```

Run ABACUS as usual, and it will automatically detect the buffer overflow problems and memory leaks. It is also possible to [use GDB with binaries built by Address Sanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizerAndDebugger).

[Valgrind](https://valgrind.org/) is another option for performing dynamic analysis.

## Generating code coverage report

This feature requires using GCC compiler. We use `gcov` and `lcov` to generate code coverage report.

1. Add `-DENABLE_COVERAGE=ON` for CMake configure command.

    ```bash
    cmake -B build -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
    ```

2. Build, install ABACUS, and run test cases. Please note that since all optimizations are disabled to gather running status line by line, the performance is drastically decreased. Set a longer time out to ensure all tests are executed.

    ```bash
    cmake --build build --target test ARGS="-V --timeout 21600"
    ```

3. Generate HTML report.

    ```bash
    cd build/
    make lcov
    ```

Now you can copy `build/lcov` to your local device, and view `build/lcov/html/all_targets/index.html`.

We use [Codecov](https://codecov.io/) to host and visualize our [**code coverage report**](https://app.codecov.io/gh/deepmodeling/abacus-develop). Analysis is scheduled after a new version releases; this [action](https://github.com/deepmodeling/abacus-develop/actions/workflows/coverage.yml) can also be manually triggered.

## Submitting a Pull Request

1. [Fork](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo) the [ABACUS repository](https://github.com/deepmodeling/abacus-develop). If you already had an existing fork, [sync](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork) the fork to keep your modification up-to-date.

2. Pull your forked repository, create a new git branch, and make your changes in it:

     ```shell
     git checkout -b my-fix-branch
     ```

3. Coding your patch, including appropriate test cases and docs.
To run a subset of unit test, use `ctest -R <test-match-pattern>` to perform tests with name matched by given pattern.

4. After tests passed, commit your changes [with a proper message](#commit-message-guidelines).

5. Push your branch to GitHub:

    ```shell
    git push origin my-fix-branch
    ```

6. In GitHub, send a pull request (PR) with `deepmodeling/abacus-develop:develop` as the base repository. It is required to document your PR following [our guidelines](#commit-message-guidelines).

7. After your pull request is merged, you can safely delete your branch and sync the changes from the main (upstream) repository:

- Delete the remote branch on GitHub either [through the GitHub web UI](https://docs.github.com/en/repositories/configuring-branches-and-merges-in-your-repository/managing-branches-in-your-repository/deleting-and-restoring-branches-in-a-pull-request#deleting-a-branch-used-for-a-pull-request) or your local shell as follows:

    ```shell
    git push origin --delete my-fix-branch
    ```

- Check out the master branch:

    ```shell
    git checkout develop -f
    ```

- Delete the local branch:

    ```shell
    git branch -D my-fix-branch
    ```

- Update your master with the latest upstream version:

    ```shell
    git pull --ff upstream develop
    ```

## Commit message guidelines

A well-formatted commit message leads a more readable history when we look through some changes, and helps us generate change log.
We follow up [The Conventional Commits specification](https://www.conventionalcommits.org) for commit message format.
This format is also required for PR title and message.
The commit message should be structured as follows:

```text
<type>[optional scope]: <description>

[optional body]

[optional footer(s)]
```

- Header
  - type: The general intention of this commit
    - `Feature`: A new feature
    - `Fix`: A bug fix
    - `Docs`: Only documentation changes
    - `Style`: Changes that do not affect the meaning of the code
    - `Refactor`: A code change that neither fixes a bug nor adds a feature
    - `Perf`: A code change that improves performance
    - `Test`: Adding missing tests or correcting existing tests
    - `Build`: Changes that affect the build system or external dependencies
    - `CI`: Changes to our CI configuration files and scripts
    - `Revert`: Reverting commits
  - scope: optional, could be the module which this commit changes; for example, `orbital`
  - description: A short summary of the code changes: tell others what you did in one sentence.
- Body: optional, providing detailed, additional, or contextual information about the code changes, e.g. the motivation of this commit, referenced materials, the coding implementation, and so on.
- Footer: optional, reference GitHub issues or PRs that this commit closes or is related to. [Use a keyword](https://docs.github.com/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue#linking-a-pull-request-to-an-issue-using-a-keyword) to close an issue, e.g. "Fix #753".

Here is an example:

```text
Fix(lcao): use correct scalapack interface.

`pzgemv_` and `pzgemm_` used `double*` for alpha and beta parameters but not `complex*` , this would cause error in GNU compiler.

Fix #753.
```
