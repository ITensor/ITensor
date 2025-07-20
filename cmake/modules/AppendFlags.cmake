macro(append_flags _flags _append_flag)

  string(STRIP "${_append_flag}" _append_flag )
  set(${_flags} "${${_flags}} ${_append_flag}")
  string(STRIP "${${_flags}}" ${_flags})

endmacro()