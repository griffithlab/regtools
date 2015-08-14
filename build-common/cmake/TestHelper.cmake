#Define integration test
macro(def_integration_test exe_tgt testName script)
    set(PPATH "${BC_PYTHONPATH_EXTRA}:${CMAKE_SOURCE_DIR}/build-common/python:$ENV{PYTHONPATH}")
    add_test(
        NAME ${testName}
        COMMAND sh -ec "PYTHONPATH='${PPATH}' ${CMAKE_CURRENT_SOURCE_DIR}/${script} $<TARGET_FILE:${exe_tgt}>"
    )
    set_tests_properties(${testName} PROPERTIES LABELS integration)
endmacro(def_integration_test testName script)
