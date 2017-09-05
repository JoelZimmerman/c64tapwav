#ifndef _TAP_H
#define _TAP_H

#define TAP_RESOLUTION 8

struct tap_header {
	char identifier[12];
	char version;
	char reserved[3];
	unsigned int data_len;
};

#endif  // !defined(_TAP_H)
