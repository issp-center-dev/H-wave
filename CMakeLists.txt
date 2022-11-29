cmake_minimum_required(VERSION 3.0.2)
project(QLMS)

enable_testing()


#foreach(test_type coulombintra coulombinter exchange hund ising pairhop pairlift)
foreach(test_type coulombintra coulombinter exchange hund pairhop pairlift)
  add_test(
    NAME UHFr_${test_type}
    COMMAND bash ${CMAKE_SOURCE_DIR}/tests/runtest.sh tests/uhfr/${test_type} ${CMAKE_SOURCE_DIR}
    )
endforeach()

#foreach(test_type CoulombIntra CoulombInter Exchange Hund Ising PairHop PairLift)
foreach(test_type CoulombIntra CoulombInter Exchange Ising PairHop PairLift)
  add_test(
    NAME UHFk_${test_type}
    COMMAND bash ${CMAKE_SOURCE_DIR}/tests/runtest.sh tests/uhfk/${test_type} ${CMAKE_SOURCE_DIR}
    )
endforeach()
