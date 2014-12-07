fun! Build()
    compiler gcc
    set makeprg=make\ -j4
    make
endf

fun! BuildCython()
    compiler cython
    set makeprg=make\ -j4\ cython
    make
endf

map <F5> :call Build()<cr>
map <s-F5> :call BuildCython()<cr>
