# Style Guide

Follow the [Pony standard library Style Guide](https://github.com/ponylang/ponyc/blob/main/STYLE_GUIDE.md) with the following expceptions.

- Use 2 empty lines to separate functions. This allows to use single lines to separate group of lines withing a function, to create logical groupings.
- Use 1 empty lines to separate attributes.
- Don't write multiple instructions on the same line. You mustn't use `;` to separate instructions.
- All functions and variables, public or not, must have a docstring.
- All public functions nust have unit tests.
- Don't return on a single line `if`/`while`/`repeat`. `return` instruction must be on a single line to see the fast exit. Like:
```pony
    if one() then
      return that
    end
```
The same applies for `error` that must be on a its own line.`
- Never remove docstring text if its content is still valid. Add complement or ask the user before the change.
- Favor using the names `i`, `j`, `k` in that order for loop variables.
