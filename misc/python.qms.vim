" .qms syntax highlighting for Neovim
" Copy to ~/.config/nvim/after/syntax/python.qms.vim
" Or symlink: ln -s $(pwd)/misc/python.qms.vim ~/.config/nvim/after/syntax/
"
" Also need ftdetect: see misc/qms-ftdetect.lua

if exists('b:current_syntax') && b:current_syntax == 'python.qms'
  finish
endif

runtime! syntax/python.vim
unlet! b:current_syntax
let b:current_syntax = 'python.qms'

" formula keyword
syntax keyword qmsFormula formula
highlight default link qmsFormula Statement

" LaTeX inside formula
syntax region qmsFormulaString start=/\%(formula\)\@<=\s\+/ end=/$/
      \ containedin=pythonStatement oneline
highlight default link qmsFormulaString String

" Section markers
syntax match qmsSectionComment /^#\s*=\+.\+=\+\s*$/
highlight default link qmsSectionComment SpecialComment
