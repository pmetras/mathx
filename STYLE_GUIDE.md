# Style Guide

bits follows the [Pony standard library Style Guide](https://github.com/ponylang/ponyc/blob/main/STYLE_GUIDE.md) with the following expceptions.

- Use 2 empty lines to separate functions. This allows to use single lines to separate group of lines withing a function, to create logical groupings.
- Use 1 empty lines to separate attributes.
- Don't write multiple instructions on the same line. We mustn't use `;` to separate instructions.
- All functions and variables, public or not, must have a docstring.
- All public functions nust have unit tests.