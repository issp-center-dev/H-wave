cmake_minimum_required(VERSION 3.0.2)
project(QLMS)

enable_testing()

add_test(
  NAME UHFr_CoulombIntra
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfr/coulombintra ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFr_CoulombInter
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfr/coulombinter ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFr_Exchange
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfr/exchange ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFr_Hund
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfr/hund ${CMAKE_SOURCE_DIR}
)
# add_test(
#   NAME UHFr_Ising
#   COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfr/ising ${CMAKE_SOURCE_DIR}
# )
add_test(
  NAME UHFr_PairHop
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfr/pairhop ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFr_PairLift
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfr/pairlift ${CMAKE_SOURCE_DIR}
)

add_test(
  NAME UHFk_CoulombIntra
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfk/CoulombIntra ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFk_CoulombInter
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfk/CoulombInter ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFk_Exchange
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfk/Exchange ${CMAKE_SOURCE_DIR}
)
# add_test(
#   NAME UHFk_Hund
#   COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfk/Hund ${CMAKE_SOURCE_DIR}
# )
add_test(
  NAME UHFk_Ising
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfk/Ising ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFk_PairHop
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfk/PairHop ${CMAKE_SOURCE_DIR}
)
add_test(
  NAME UHFk_PairLift
  COMMAND bash ${CMAKE_SOURCE_DIR}/test/runtest.sh test/uhfk/PairLift ${CMAKE_SOURCE_DIR}
)

