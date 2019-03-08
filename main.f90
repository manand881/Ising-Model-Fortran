program ising !2D Monte Carlo Simultion Of Ising Model
    !Lisa Larrimore,lisal@succs.swarthmore.edu
    !3 May 2002
    !Physics 114 Final Project

    !This program is adapted from the Ising Model program written in
    !BASIC by Elaine Chandler that appears on p.184 of David Chandler's
    !Introduction to Modern Statistical Mechanics.

    !The Input parameters of this Program are in "ising.in", and they
    !allow the size, length and initial configuration of the simulation
    !to be changed. See comments in file.

    !This Program has three output files:
    !
    !"Sping-array"      Contains snapshots of the spin lattice at the end of
    !                   each temperature run ( or throughout the middle of the
    !                   run, if only looking at one temperature). Can be
    !                   visualized with the IDL program see_spins.pro
    !
    !"magnetization"    Contains four columns: each temperature,the
    !                   average magnetization at that temp,the avg magnetization
    !                   squared at that temp, and the susceptibility.
    !
    !"energy"           contains four columns: each temperature, the
    !                   average energy at that temp,the avg energy squared
    !                   at that temp,and the heat capacity.
    !

    implicit none

    
end program ising