/*
 * NSDefinitions.h
 *
 */
#ifndef NSDEFINITIONS_H_
#define NSDEFINITIONS_H_

/* boundary values */
#define NO_SLIP		1
#define FREE_SLIP	2
#define OUTFLOW 	3

/**
 * Boundary cells
 */
#define P_L 32  /* 0b0100000*/
#define P_R 64  /* 0b1000000*/

/* cell values fluid/obstacle (Bits: center|east|west|south|north),
 * 1 if corresponding cell is a fluid cell, 0 if it is an obstacle */
#define C_F			255
#define C_B			1

/**
 * Boundary cell classification |center|east|west|south|north| (1 = fluid)
 */
#define B_N			1  /* 0b0000001*/
#define B_S			2  /* 0b0000010*/
#define B_W			4  /* 0b0000100*/
#define B_O			8  /* 0b0001000*/
#define B_C 		16 /* 0b0010000*/
#define B_NW		5  /* 0b0000101*/
#define B_SW		6  /* 0b0000110*/
#define B_NO		9  /* 0b0001001*/
#define B_SO		10 /* 0b0001010*/

#endif
