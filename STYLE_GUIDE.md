# Style Guide

**Important**: After every major change to the code, you must check that the tests pass and that the documentation is up-to-date. Enventually, update existing tests or add new ones to have a complete coverage.

Follow the [Pony standard library Style Guide](https://github.com/ponylang/ponyc/blob/main/STYLE_GUIDE.md) with the following expceptions.

- Use 2 empty lines to separate functions. This allows to use single lines to separate group of lines withing a function, to create logical groupings.
- Use 1 empty lines to separate attributes.
- Don't write multiple instructions on the same line. You mustn't use `;` to separate instructions.
- All classes, agents, primitives, types, functions, fields and variables, etc., public or not, must have a docstring. Everything that supports the docstring syntax should have one documenting its role.
- Docstrings must be written in good English, without shortcuts. They must be understandable by any developer, not experts. Java or Python documentations are good examples of well written documentation.
- When a docstring includes a reference to another method of variable documentation (i.e. See other_fun), this reference must be a link to the other documentation.
- All public functions nust have unit tests.
- `return` instruction must be on a single line to see the fast exit. The same applies for `error` that must be on a its own line. Like:
```pony
    if one() then
      return that
    else
      error
    end
```
- Never remove docstring text if its content is still valid. Add complement or ask the user before changing it.
- Use `Range` instead of `while` for simple loops. The compiler will optimize it and it's simpler to read. Also, the scope of the loop variable is limited to the loop.
- Favor using the names `i`, `j`, `k` in that order for loop variables when they are general. If they have a meaning, try to use a name that carries that meaning. For example, if the variable is an item counter and that the item is central to the algorithm, use `item_counter` instead.
- When an `if` expression result is assigned to a variable, the `let` (if present), `if` and `then` keyword should be on the same line, like
```pony
      term = if d1 <= 255 then
          (let q, _) = term._short_div(d1.u8())
          q
        else
          _div(term, MPFRep.from_f64(d1.f64(), p2), p2)
        end
```
- `if` instructions must have `if` and `then` in the same line if possible, `else` on another line if present, and `end` on another line.
- All instructions must be on their own line if possible. Don't use `;` to put 2 instrunctions on the same line.
- `match` cases must be on different lines.
