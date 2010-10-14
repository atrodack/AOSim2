/*** shutters.c ************************************************************
 * 
 * This is the command-line interface for the shutters and lasers.
 * 
 * The Linux version is currently built on top of the PCI-DIO24 device driver 
 * by Warren J. Jasper at North Carolina State University (wjasper@tx.ncsu.edu).  
 * That driver is GPL, and therefore so is this.  
 * The latest source for the PCI-DIO24 driver is located at 
 * ftp://lx10.tx.ncsu.edu/pub/Linux/drivers.
 * 
 * ***********************************************************************
 * Written by Johanan L. Codona
 * Center for Astronomical Adaptive Optics
 * Steward Observatory
 * University of Arizona
 * Tucson, AZ 85721
 * jcodona@as.arizona.edu
 * ***********************************************************************
 * CHANGES:
 * 20061221: JLCodona: Ver 0.1.  Based this on my old uDM140 driver.  
 *           Not a shared library this time, just a simple command.
 * ***************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/ioctl.h>

#include <unistd.h>
#include <string.h>
#include <stdint.h>
#include <sys/ioctl.h>
#include <getopt.h>
		
#include <pci-dio24.h>
		
#define DEV_A "/dev/daq/DIO24/A"
#define DEV_B "/dev/daq/DIO24/B"
#define DEV_C "/dev/daq/DIO24/C"

/* Return codes */
#define NO_ERROR         0
#define DEV_ERROR       -1
#define IO_ERROR        -2
#define PROGRAM_ERROR   -3

typedef int boolean;
#define false (0)
#define true  (-1)

#define DEBUG false
#define STIME 500000

#define COMMAND_OPTIONS "?crxOIS:01L:ocA:B:C:Z:D:d:"

static int A=-1, B=-1, C=-1;

static int debug = 0;

int wrA(int x)
{
	int retval;
	
	if(debug) fprintf(stderr,"A <-- 0x%04X (%d)\n", x,x);
	retval = write(A,&x,1);
	return retval;
}

int wrB(int x)
{
	int retval;
	
	if(debug) fprintf(stderr,"B <-- 0x%04X (%d)\n", x,x);
	retval = write(B,&x,1);
	return retval;
}

int wrC(int x)
{
	int retval;
	
	if(debug) fprintf(stderr,"C <-- 0x%04X (%d)\n", x,x);
	retval = write(C,&x,1);
	return retval;
}



/*** open_dio24 **********************************************************
 * Open the device nodes for operation and package handles.
 **********************************************************************/
int open_dio24(char *nodeA, char *nodeB, char *nodeC)
{
	int err; 
	char buf[120];

	// fprintf(stderr,"dm_open(%s,%s,%s).\n",nodeA,nodeB,nodeC);
	
	// open the device nodes
  	if((A = open(nodeA, O_RDWR))==-1) {
    		sprintf(buf,"error opening node A: %s", nodeA);
    		perror(buf);
		return DEV_ERROR;
  	}
 
	if((B = open(nodeB, O_RDWR))==-1) {
    		sprintf(buf,"error opening node B: %s", nodeB);
		perror(buf);
		return DEV_ERROR;
	}
 
  	if((C = open(nodeC, O_RDWR))==-1) {
 		sprintf(buf,"error opening node C: %s", nodeC);
   		perror(buf);
		return DEV_ERROR;
	}
	
	return NO_ERROR;
}

/*** dm_close *********************************************************
 * Close the device nodes for operation and package handles.
 **********************************************************************/
int close_dio24(void)
{
	fprintf(stderr,"dm_close\n");

	if(A >= 0) close(A);
	if(B >= 0) close(B);
	if(C >= 0) close(C);

	return NO_ERROR;
}

/*** configIO *******************************************************
 * Set for output, etc.
 **********************************************************************/

int configOUTPUTS(void)
{
	int err;
		
	fprintf(stderr,"configOUTPUTS\n");
	
	// set the registers to be outputs....
	if(ioctl(A, DIO_SET_DIRECTION, PORT_OUTPUT)==-1) {
		perror("error configuring node A");
		return DEV_ERROR;
	}
	
	if(ioctl(B, DIO_SET_DIRECTION, PORT_OUTPUT)==-1) {
		perror("error configuring node B");
			return DEV_ERROR;
		}
		
// Set this to be outputs too when you need it. 
		
// 	if(ioctl(DMI.C, DIO_SET_DIRECTION, LOW_PORT_OUTPUT)==-1) {
// 		perror("error configuring node CL");
// 			return DEV_ERROR;
// 		}
// 		
// 	if(ioctl(DMI.C, DIO_SET_DIRECTION, HIGH_PORT_OUTPUT)==-1) {
// 		perror("error configuring node CH");
// 		return DEV_ERROR;
// 		}
		
	return NO_ERROR;
}

int configINPUTS(void)
{
	int err;
		
	fprintf(stderr,"configINPUTS (This is not the normal configuration.)\n");
	
	// set the registers to be outputs....
	if(ioctl(A, DIO_SET_DIRECTION, PORT_INPUT)==-1) {
		perror("error configuring node A");
		return DEV_ERROR;
	}
	
	if(ioctl(B, DIO_SET_DIRECTION, PORT_INPUT)==-1) {
		perror("error configuring node B");
		return DEV_ERROR;
	}
		
// Set this to be outputs too when you need it. 
		
// 	if(ioctl(DMI.C, DIO_SET_DIRECTION, LOW_PORT_INPUT)==-1) {
// 		perror("error configuring node CL");
// 			return DEV_ERROR;
// 		}
	// 		
// 	if(ioctl(DMI.C, DIO_SET_DIRECTION, HIGH_PORT_INPUT)==-1) {
// 		perror("error configuring node CH");
// 		return DEV_ERROR;
// 		}
		
	return NO_ERROR;
}


int main(int argc, char **argv)
{
	int LASER=-1;
	int SHUTTER=-1;
	
	int value;
	
	int c;
	int digit_optind = 0;
	
	// open device
	int err;
	if((err = open_dio24(DEV_A, DEV_B, DEV_C)) != NO_ERROR) {
		fprintf(stderr,"Error opening PCI-DIO24 device nodes: err = %d (0x%x)\n",err,err);
		exit -1;
	}
	
	// parse args and execute commands in order
	while (1) {
		int this_option_optind = optind ? optind : 1;
		int option_index = 0;
		static struct option long_options[] = {
			{"config", no_argument, 0, 0},
			{"dark",   no_argument, 0, 0},
			{"ref",    no_argument, 0, 0},
			{"no_ref", no_argument, 0, 0},
			{"halo",   no_argument, 0, 0},
			{"no_halo",no_argument, 0, 0},
			{"red",    no_argument, 0, 0},
			{"green",  no_argument, 0, 0},
			{"yellow", no_argument, 0, 0},
			{"inputs", no_argument, 0, 0},
			{"outputs", no_argument, 0, 0},
			{0, 0, 0, 0}
		};

		c = getopt_long(argc, argv, COMMAND_OPTIONS, long_options, &option_index);
		if (c == -1)
			break;

		switch (c) {
			case 0:
				{
				char *larg = long_options[option_index].name;

				/* if (debug) {
					printf ("found long option: %s", larg);
					if (optarg)
						printf (" with arg %s", optarg);
					printf ("\n");
				}*/
				
				// LASER COMMANDS
				if(!strcmp(larg,"red")) {
					if(debug) printf("Turning RED laser ON.\n");
					wrA(8);

				} else if(!strcmp(larg,"green")) {
					if(debug) printf("Turning GREEN laser ON.\n");
					wrA(2);

				} else if(!strcmp(larg,"yellow")) {
					if(debug) printf("Turning YELLOW laser ON.\n");
					wrA(1);

				} else if(!strcmp(larg,"dark")) {
					if(debug) printf("Turning ALL lasers OFF.\n");
					wrA(0);

				// REFERENCE BEAM COMMANDS
				} else if(!strcmp(larg,"ref")) {
					int shutter=0;
					int START=2<<(2*shutter);

					if(debug) printf("Turning REFERENCE beam ON.\n");
					wrB(START); usleep(STIME); wrB(0);
					if(debug) printf("Done.\n");

				} else if(!strcmp(larg,"no_ref")) {
					int shutter=0;
					int START=1<<(2*shutter);

					if(debug) printf("Turning REFERENCE beam ON.\n");
					wrB(START); usleep(STIME); wrB(0);
					if(debug) printf("Done.\n");

				// HALO PATH COMMANDS
				} else if(!strcmp(larg,"halo")) {
					int shutter=1;
					int START=2<<(2*shutter);

					if(debug) printf("Turning HALO path ON.\n");
					wrB(START); usleep(STIME); wrB(0);
					if(debug) printf("Done.\n");

				} else if(!strcmp(larg,"no_halo")) {
					int shutter=1;
					int START=1<<(2*shutter);

					if(debug) printf("Turning HALO path OFF.\n");
					wrB(START); usleep(STIME); wrB(0);
					if(debug) printf("Done.\n");

				// DIO COMMANDS
				} else if(!strcmp(larg,"outputs")) {
					if(debug) printf ("configuring interface for OUTPUT.\n");
					configOUTPUTS();

				} else if(!strcmp(larg,"inputs")) {
					if(debug) printf ("Resetting the interface for INPUTS (NOTE BENE!).\n");
					configINPUTS();

				} 

				
				}
				break;
				
			case 'O':
				if(debug) printf ("configuring interface for OUTPUT.\n");
				configOUTPUTS();
				break;

			case 'r':
			case 'I':
			case 'x':
				if(debug) printf ("Resetting the interface for INPUTS (NOTE BENE!).\n");
				configINPUTS();
				break;
				
			case 'b':
				if(debug) printf ("option b\n");
				break;

			case 'z':
			case 'Z':
				if(debug) printf ("Sleeping for %d usecs.\n",atoi(optarg));
				usleep(atoi(optarg));
				break;

			case 'S':
				SHUTTER = atoi(optarg);
				if(debug) printf ("SHUTTER %d\n", SHUTTER);
				break;

			case 'L':
				LASER = atoi(optarg);
				if(debug) printf ("LASER %d\n", LASER);
				break;

			case 'A':
				value = atoi(optarg);
				if(debug) printf("Register A <- %d\n", value);
				wrA(value);
				break;

			case 'd':
			case 'D':
				debug = atoi(optarg);
				break;

			case 'B':
				value = atoi(optarg);
				if(debug) printf("Register B <- %d\n", value);
				wrB(value);
				break;

			case 'C':
				value = atoi(optarg);
				if(debug) printf("Register C <- %d\n", value);
				wrC(value);
				break;

			case 'c':
				if(SHUTTER<0) 
					fprintf(stderr,"ERROR: Close WHICH shutter?\n");
				else {
					int START=1<<(2*SHUTTER);
					int STOP=0;
					
					wrB(START);
					if(debug) printf("Closing Shutter %d... ", SHUTTER);
					usleep(STIME);
					wrB(STOP);
					if(debug) printf("Done.\n");
				}

				break;
				
			case 'o':
				if(SHUTTER<0) 
					fprintf(stderr,"ERROR: Open WHICH shutter?\n");
				else {
					int START=2<<(2*SHUTTER);
					int STOP=0;
					
					wrB(START);
					if(debug) printf("Opening Shutter %d (0x%x)... ",SHUTTER,START);
					usleep(STIME);
					wrB(STOP);
					if(debug) printf("Done.\n");
				}

				break;
				
			case '1':
				if(LASER<0) 
					fprintf(stderr,"ERROR: Turn WHICH laser on?\n");
				else {
					wrA(1<<LASER);
					if(debug) printf("Activating laser %d... ", LASER);
				}
				break;
				
			case '0':
				wrA(0);
				if(debug) 
					printf("Turning ALL A-register devices off.\n");
				break;
				
			case '?':
				fprintf(stderr, "options: %s\n", COMMAND_OPTIONS);
				break;

			default:
				fprintf (stderr,"?? getopt returned character code 0%o ??\n", c);
		}
	}

	if (optind < argc) {
		fprintf (stderr, "non-option ARGV-elements: ");
		while (optind < argc)
			fprintf (stderr,"%s ", argv[optind++]);
		fprintf (stderr,"\n");
	}

	exit (0);
}



