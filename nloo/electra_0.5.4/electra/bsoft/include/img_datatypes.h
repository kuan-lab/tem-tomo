/*
	img_datatypes.h
	Functions to check and convert image datatypes
	Author: Bernard Heymann
	Created: 19990321 	Modified: 20040504
*/

#include "rwimg.h"

// Function prototypes
unsigned long   gettypesize(DataType type);
DataType	getdatatype(char letter);
DataType	get_datatype_from_string(char* string);
char* 		get_string_from_datatype(DataType type);
char* 		get_string_from_values(char* data, int bit, int channels, DataType type);
double 		get_value_from_datatype(char* data, DataType type);
unsigned char*  set_value_datatype(float value, DataType type, int channels);
int			img_to_datatype(Bimage* p, DataType newtype);
int			img_to_bit(Bimage* p);
int			img_bit_to_byte(Bimage* p);
int			img_to_byte(Bimage* p);
int			img_to_signed_byte(Bimage* p);
int			img_to_unsigned_short(Bimage* p);
int			img_to_short(Bimage* p);
int			img_to_int(Bimage* p);
int			img_to_float(Bimage* p);
int 		img_gray2RGB(Bimage* p);
int 		img_gray2RGBA(Bimage* p);
int 		img_gray2color(Bimage* p, int alpha);
int 		img_RGB2gray(Bimage* p);
int 		img_indexed2gray(Bimage* p);
int 		img_indexed2RGB(Bimage* p);
int			img_add_alpha(Bimage* p);
int			img_remove_alpha(Bimage* p);

