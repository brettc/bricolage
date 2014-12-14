fun! Build()
    compiler gcc
    set makeprg=make\ -j4
    make
endf

fun! BuildCython()
    compiler cython
    echo "set makeprg=make ".expand("%:p:r").cpp
    " exe "set makeprg=make ".expand("%:p:r").cpp
    make
endf

map <F10> :call Build()<cr>
map <s-F10> :call BuildCython()<cr>
