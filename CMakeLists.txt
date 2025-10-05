cmake_minimum_required(VERSION 3.30)
project(TFM_UPC)

set(CMAKE_CXX_STANDARD 20)

add_executable(Esercizio1 Heat_Conduction_1D.cpp
        Functions/Functions.cpp)

add_executable(Esercizio2 Heat_Conduction_2D.cpp
        Functions/Functions.cpp)

add_executable(Esercizio3 Smith_Hutton.cpp)

add_executable(Esercizio3.2 Smith_Hutton_Upgraded.cpp)

add_executable(Esercizio4 Lid_Driven_Cavity.cpp)

add_executable(Esercizio5 Differentially_Heated_Cavity.cpp)

add_executable(Esercizio6 DHC_Turbulent.cpp)