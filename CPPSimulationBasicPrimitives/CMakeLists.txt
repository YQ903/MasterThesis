cmake_minimum_required(VERSION 3.6)

project(GeometricPrimitives)
include_directories(SYSTEM /usr/include/eigen3)

set(DQ_PREFIX /usr/local)
include_directories(SYSTEM ${DQ_PREFIX}/include)
find_package (Threads REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package (osqp REQUIRED)
find_package (OsqpEigen REQUIRED)
#find_package(catkin REQUIRED COMPONENTS
#  OsqpEigen
#)




add_executable(inversekinematic
    src/InverseKinematic.cpp
    src/FrankaRobot.cpp
    src/RealSenseRobot.cpp
    src/geomJac.cpp
    src/convert_matrix2csv.cpp
    src/convert_csv2matrix.cpp
    include/convert_matrix2csv.h
    include/convert_csv2matrix.h
    include/RealSenseRobot.h
    include/geomJac.h
    include/FrankaRobot.h)

TARGET_LINK_LIBRARIES(inversekinematic
    pthread
    ${DQ_PREFIX}/lib/libdqrobotics.so
    ${DQ_PREFIX}/lib/libdqrobotics-interface-vrep.so
    yaml-cpp
)



add_executable(quadraticprogramming
    src/QuadraticProgramming.cpp
    src/FrankaRobot.cpp
    src/RealSenseRobot.cpp
    src/geomJac.cpp
    src/convert_matrix2csv.cpp
    src/convert_csv2matrix.cpp
    include/convert_matrix2csv.h
    include/convert_csv2matrix.h
    include/RealSenseRobot.h
    include/geomJac.h
    include/FrankaRobot.h)



TARGET_LINK_LIBRARIES(quadraticprogramming 
    pthread
    dqrobotics
    dqrobotics-interface-vrep
    osqp::osqp
    OsqpEigen::OsqpEigen
    yaml-cpp
)


add_executable(maniqphuman
    src/ManiQPHuman.cpp
    src/FrankaRobot.cpp
    src/RealSenseRobot.cpp
    src/geomJac.cpp
    src/convert_matrix2csv.cpp
    src/convert_csv2matrix.cpp
    include/convert_matrix2csv.h
    include/convert_csv2matrix.h
    include/RealSenseRobot.h
    include/geomJac.h
    include/FrankaRobot.h
    include/logmap.h
    include/symm2vec.h
    include/spdToVec.h
    include/jacobianEst.h
    include/tmprod.h
    include/readData.h
    include/manipulabilityJacobian.h
    include/jacobianEstVector.h
    
    src/convert_matrix2csv.cpp
    src/convert_csv2matrix.cpp
    src/RealSenseRobot.cpp
    src/geomJac.cpp
    src/FrankaRobot.cpp
    src/logmap.cpp
    src/symm2vec.cpp
    src/spdToVec.cpp
    src/jacobianEst.cpp
    src/tmprod.cpp
    src/readData.cpp
    src/manipulabilityJacobian.cpp
    src/jacobianEstVector.cpp
    src/redManipulabilityJacobian.cpp

    )



TARGET_LINK_LIBRARIES(maniqphuman
    pthread
    dqrobotics
    dqrobotics-interface-vrep
    osqp::osqp
    OsqpEigen::OsqpEigen
    yaml-cpp
)

add_executable(maniqp
    include/franka_ik_He.hpp
    src/fik.cpp
    include/fik.h

    src/ManiQP.cpp
    src/FrankaRobot.cpp
    src/RealSenseRobot.cpp
    src/geomJac.cpp
    src/convert_matrix2csv.cpp
    src/convert_csv2matrix.cpp
    include/convert_matrix2csv.h
    include/convert_csv2matrix.h
    include/RealSenseRobot.h
    include/geomJac.h
    include/FrankaRobot.h
    include/logmap.h
    include/symm2vec.h
    include/spdToVec.h
    include/jacobianEst.h
    include/tmprod.h
    include/readData.h
    include/manipulabilityJacobian.h
    include/jacobianEstVector.h
    
    src/convert_matrix2csv.cpp
    src/convert_csv2matrix.cpp
    src/RealSenseRobot.cpp
    src/geomJac.cpp
    src/FrankaRobot.cpp
    src/logmap.cpp
    src/symm2vec.cpp
    src/spdToVec.cpp
    src/jacobianEst.cpp
    src/tmprod.cpp
    src/readData.cpp
    src/manipulabilityJacobian.cpp
    src/jacobianEstVector.cpp
    src/redManipulabilityJacobian.cpp

    )



TARGET_LINK_LIBRARIES(maniqp
    pthread
    dqrobotics
    dqrobotics-interface-vrep
    osqp::osqp
    OsqpEigen::OsqpEigen
    yaml-cpp
)
