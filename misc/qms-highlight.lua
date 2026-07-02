# .qms Syntax Highlighting for Neovim
#
# Quick install:
#   cp misc/qms-highlight.lua ~/.config/nvim/plugin/qms.lua
#
# What it does:
#   1. Treats .qms files as Python (treesitter + LSP work automatically)
#   2. Highlights the 'formula' keyword in Statement color

-- ~/.config/nvim/plugin/qms.lua
local group = vim.api.nvim_create_augroup("qms_ft", { clear = true })

-- 1. Filetype detection
vim.filetype.add({ extension = { qms = "python" } })

-- 2. Highlight 'formula' keyword in .qms files
vim.api.nvim_create_autocmd({ "BufRead", "BufNewFile", "FileType" }, {
  pattern = "*.qms",
  group = group,
  callback = function()
    vim.fn.matchadd("Statement", "\\v<formula>")
  end,
})
