all
# all rules enabled with defaults until disabled explicitly.
# all options are set explicitly, even if the value is the default value

# MD001 - Header levels should only increment by one level at a time
# MD002 - First header should be a top level header
rule 'MD002', :level => 1
# MD003 - Header style
rule 'MD003', :style => :consistent
# MD004 - Unordered list style
rule 'MD004', :style => :sublist
# MD005 - Inconsistent indentation for list items at the same level
# MD006 - Consider starting bulleted lists at the beginning of the line
# MD007 - Unordered list indentation
rule 'MD007', :indent => 3
# MD009 - Trailing spaces
exclude_rule 'MD009'
# MD010 - Hard tabs
rule 'MD010', :exclude_code_blocks => true
# MD011 - Reversed link syntax
# MD012 - Multiple consecutive blank lines
exclude_rule 'MD012'
# MD013 - Line length
exclude_rule 'MD013'
# MD014 - Dollar signs used before commands without showing output
exclude_rule 'MD014'
# MD018 - No space after hash on atx style header
# MD019 - Multiple spaces after hash on atx style header
# MD020 - No space inside hashes on closed atx style header
# MD021 - Multiple spaces inside hashes on closed atx style header
# MD022 - Headers should be surrounded by blank lines
# MD023 - Headers must start at the beginning of the line
# MD024 - Multiple headers with the same content
rule 'MD024', :allow_different_nesting => true
# MD025 - Multiple top level headers in the same document
rule 'MD025', :level => 1
# MD026 - Trailing punctuation in header
exclude_rule 'MD026' # , :punctuation => ".,;:!?"
# MD027 - Multiple spaces after blockquote symbol
# MD028 - Blank line inside blockquote
# MD029 - Ordered list item prefix
rule 'MD029', :style => :ordered
# MD030 - Spaces after list markers
rule 'MD030', :ul_single => 1, :ol_single => 1, :ul_multi => 1, :ol_multi => 1
# MD031 - Fenced code blocks should be surrounded by blank lines
# MD032 - Lists should be surrounded by blank lines
# MD033 - Inline HTML
# MD034 - Bare URL used
# MD035 - Horizontal rule style
rule 'MD035', :style => :consistent
# MD036 - Emphasis used instead of a header
rule 'MD036', :punctuation => ".,;:!?"
# MD037 - Spaces inside emphasis markers
# MD038 - Spaces inside code span elements
# MD039 - Spaces inside link text
# MD040 - Fenced code blocks should have a language specified
exclude_rule 'MD040'
# MD041 - First line in file should be a top level header
exclude_rule 'MD041' #, :level => 1
# MD046 - Code block style
rule 'MD046', :style => :fenced
# MD047 - File should end with a single newline character

