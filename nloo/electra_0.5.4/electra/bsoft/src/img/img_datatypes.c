/*
	img_datatypes.c
	Functions to check and convert image datatypes
	Author: Bernard Heymann
	Created: 19990321 	Modified: 20050101
*/

#include "rwimg.h"
#include "img_datatypes.h"
#include "img_symmetry.h"
#include "img_complex.h"
#include "utilities.h"

// Declaration of global variables
extern int 	verbose;		// Level of output to the screen
extern long memory;			// Total memory allocated 

/************************************************************************
@Function: gettypesize
@Description:
	Get the size of a datatype.
@Algorithm:
	This function is used for calculating image data sizes for allocating
	memory and reading the data.
@Arguments:
	DataType type	data type (defined in rwimg.h).
@Returns:
	unsigned long   size of data type, if < 0 the data type is not supported.
**************************************************************************/
unsigned long   gettypesize(DataType type)
{
	unsigned long   size;
	
	switch ( type ) {
		case UChar: case SChar:  size = sizeof(char); break;
		case UShort: case Short: size = sizeof(short); break;
		case Int: 			 	 size = sizeof(int); break;
		case Float: 			 size = sizeof(float); break;
		case Double: 			 size = sizeof(double); break;
		case ComplexShort:  	 size = sizeof(complex_short); break;
		case ComplexInt:  	 	 size = sizeof(complex_int); break;
		case ComplexFloat:  	 size = sizeof(complex_float); break;
		case ComplexDouble:  	 size = sizeof(complex_double); break;
		case Polar:  	 		 size = sizeof(polar); break;
		default: size = 0;
	}
	
	return(size);
}

/************************************************************************
@Function: getdatatype
@Description:
	Get the data type indicated by a single letter code.
@Algorithm:
	This function is used in optional command-line arguments to indicate 
	a new data type for an image.
@Arguments:
	char letter 	letter indicating data type.
@Returns:
	DataType type	data type (defined in rwimg.h).
**************************************************************************/
DataType	getdatatype(char letter)
{
	DataType		type;
	
	switch ( letter ) {
		case 'b': type = UChar; break;
		case 'c': type = SChar; break;
		case 'u': type = UShort; break;
		case 's': type = Short; break;
		case 'i': type = Int; break;
		case 'f': type = Float; break;
		case 'd': type = Double; break;
		case 'S': type = ComplexShort; break;
		case 'I': type = ComplexInt; break;
		case 'F': type = ComplexFloat; break;
		case 'D': type = ComplexDouble; break;
		case 'P': type = Polar; break;
		default: type = Unknown_Type; break;
	}
	
	return(type);
}

/************************************************************************
@Function: get_datatype_from_string
@Description:
	Get the data type from a string.
@Algorithm:
	.
@Arguments:
	char* string 	string describing the data type.
@Returns:
	DataType type	data type (defined in rwimg.h).
**************************************************************************/
DataType	get_datatype_from_string(char* string)
{
	unsigned long   i;
	DataType		type = Unknown_Type;
	
	for ( i=0; i<strlen(string); i++ ) string[i] = tolower(string[i]);
	
	if ( strlen(string) == 1 ) {
		type = getdatatype(string[0]);
	} else {
		if ( strstr(string, "complex") ) {
			if ( strstr(string, "short") )
				type = ComplexShort;
			else if ( strstr(string, "int") )
				type = ComplexInt;
			else if ( strstr(string, "float") )
				type = ComplexFloat;
			else if ( strstr(string, "double") )
				type = ComplexDouble;
		} else {
			if ( strstr(string, "unsigned") ) {
				if ( strstr(string, "char") )
					type = UChar;
				else if ( strstr(string, "short") )
					type = UShort;
			} else {
				if ( strstr(string, "byte") )
					type = UChar;
				else if ( strstr(string, "char") )
					type = SChar;
				else if ( strstr(string, "short") )
					type = Short;
				else if ( strstr(string, "int") )
					type = Int;
				else if ( strstr(string, "float") )
					type = Float;
				else if ( strstr(string, "double") )
					type = Double;
				else if ( strstr(string, "polar") )
					type = Polar;
			}
		}
	}
	
	return(type);
}

/************************************************************************
@Function: get_string_from_datatype
@Description:
	Get the string representation of a datatype.
@Algorithm:
	.
@Arguments:
	DataType type	data type (defined in rwimg.h).
@Returns:
	char*			string representation of a datatype.
**************************************************************************/
char* 		get_string_from_datatype(DataType type)
{
	char*		temp = (char *) balloc(20);
	
	switch ( type ) {
		case UChar:             strcpy(temp, "unsigned char"); break;
		case SChar:             strcpy(temp, "signed char"); break;
		case UShort:            strcpy(temp, "unsigned short"); break;
		case Short:             strcpy(temp, "short"); break;
		case Int: 			 	strcpy(temp, "int"); break;
		case Float: 			strcpy(temp, "float"); break;
		case Double: 			strcpy(temp, "double"); break;
		case ComplexShort:  	strcpy(temp, "complex short"); break;
		case ComplexInt:  	 	strcpy(temp, "complex int"); break;
		case ComplexFloat:  	strcpy(temp, "complex float"); break;
		case ComplexDouble:  	strcpy(temp, "complex double"); break;
		case Polar:  	 		strcpy(temp, "polar"); break;
		default: strcpy(temp, "unknown");
	}
	
	char*		string = (char *) balloc(strlen(temp) + 1);
	strcpy(string, temp);
	bfree(temp, 20);
	
	return(string);
}

/************************************************************************
@Function: get_string_from_values
@Description:
	Get a formatted string from a pointer to data of a given type.
@Algorithm:
	.
@Arguments:
	char* data	 	the data pointer.
	int bit			bit offset if a bitmap (ignored if < 0).
	int channels	number of channels.
	DataType type	data type (defined in rwimg.h).
@Returns:
	char*			the formatted string.
**************************************************************************/
char* 		get_string_from_values(char* data, int bit, int channels, DataType type)
{
	char		astring[128] = "";
	char*		aptr = astring;
	int			c;
	
	if ( bit >  -1 ) {
		if ( 0x80 & ( data[0] << bit ) ) {
			sprintf(aptr, "1 ");
		} else {
			sprintf(aptr, "0 ");
		}
	} else for ( c=0; c<channels; c++ ) {
		switch ( type ) {
			case UChar:
			case SChar:
			case UShort:
			case Short:
			case Int:
			case Long:
				sprintf(aptr, "%.0f ", get_value_from_datatype(data, type));
				break;
			case Float:
			case Double:
				sprintf(aptr, "%g ", get_value_from_datatype(data, type));
				break;
			case ComplexShort:
				sprintf(aptr, "%.0f %.0f ", get_value_from_datatype(data, type),
					get_value_from_datatype(data+sizeof(short), type));
				break;
			case ComplexInt:
				sprintf(aptr, "%.0f %.0f ", get_value_from_datatype(data, type),
					get_value_from_datatype(data+sizeof(int), type));
				break;
			case ComplexFloat:
				sprintf(aptr, "%g %g ", get_value_from_datatype(data, type),
					get_value_from_datatype(data+sizeof(float), type));
				break;
			case ComplexDouble:
				sprintf(aptr, "%g %g ", get_value_from_datatype(data, type),
					get_value_from_datatype(data+sizeof(double), type));
				break;
			case Polar:
				sprintf(aptr, "%g %g ", get_value_from_datatype(data, type),
					get_value_from_datatype(data+sizeof(float), type));
				break;
		}
		data += gettypesize(type);
		aptr = astring + strlen(astring);
	}
	
	if ( strlen(astring) < 1 ) return(NULL);
	
	char*		thestring = (char *) balloc(strlen(astring)+1);
	strncpy(thestring, astring, strlen(astring));
	
	return(thestring);
}

/************************************************************************
@Function: get_value_from_datatype
@Description:
	Get a floating point value from a pointer to data of a given type.
@Algorithm:
	.
@Arguments:
	char* data	 	the data pointer.
	DataType type	data type (defined in rwimg.h).
@Returns:
	double			the floating point value.
**************************************************************************/
double 		get_value_from_datatype(char* data, DataType type)
{
	float 			value = 0;
	
	complex_short	csval;
	complex_int 	cival;
	complex_float 	cfval;
	complex_double 	cdval;
	
	switch ( type ) {
		case UChar:  value = *((unsigned char *) data); break;
		case SChar:  value = *((signed char *) data); break;
		case UShort: value = *((unsigned short *) data); break;
		case Short:  value = *((short *) data); break;
		case Int:    value = *((int *) data); break;
		case Float:  value = *((float *) data); break;
		case Double: value = *((double *) data); break;
		case ComplexShort:
			csval = *((complex_short *) data);
			value = sqrt(1.0*csval.re*csval.re + csval.im*csval.im);
			break;
		case ComplexInt:
			cival = *((complex_int *) data);
			value = sqrt(1.0*cival.re*cival.re + cival.im*cival.im);
			break;
		case ComplexFloat:
			cfval = *((complex_float *) data);
			value = sqrt(cfval.re*cfval.re + cfval.im*cfval.im);
			break;
		case ComplexDouble:
			cdval = *((complex_double *) data);
			value = sqrt(cdval.re*cdval.re + cdval.im*cdval.im);
			break;
		case Polar:  value = *((float *) data);  break;
		default: ;
	}
	
	return(value);
}

/************************************************************************
@Function: set_value_datatype
@Description:
	Set a pointer to a value with a given data type.
@Algorithm:
	A pointer to a memory location with the given data type is returned.
	The size of the allocated memory is channels*typesize.
@Arguments:
	float value 	the value.
	DataType type	data type (defined in rwimg.h).
	int channels	number of channels.
@Returns:
	char*			pointer to the value in the desired data type.
**************************************************************************/
unsigned char*  set_value_datatype(float value, DataType type, int channels)
{
	if ( channels < 1 ) channels = 1;
	
	int		size = (int) (channels*gettypesize(type));
	if ( size < 0 ) return(NULL);
	
	int 			i;
	unsigned char*  val_ptr = (unsigned char *) balloc(size);
	
	signed char		cval = (signed char)    value;
	unsigned short	uval = (unsigned short) value;
	short			sval = (short)          value;
	int 			ival = (int)            value;
	float 			fval = (float)          value;
	double 			dval = (double)         value;
	complex_short	csval = {(short)  value, (short)  value};
	complex_int 	cival = {(int)    value, (int)    value};
	complex_float 	cfval = {(float)  value, (float)  value};
	complex_double 	cdval = {(double) value, (double) value};
	polar		 	pval =  {(float)  value, 0};
	
	switch ( type ) {
		case UChar:
			for ( i=0; i<channels; i++ )
				val_ptr[i] = (unsigned char) value;
			break;
		case SChar:  memcpy(val_ptr, &cval, size); break;
		case UShort: memcpy(val_ptr, &uval, size); break;
		case Short:  memcpy(val_ptr, &sval, size); break;
		case Int:    memcpy(val_ptr, &ival, size); break;
		case Float:  memcpy(val_ptr, &fval, size); break;
		case Double: memcpy(val_ptr, &dval, size); break;
		case ComplexShort:  memcpy(val_ptr, &csval, size); break;
		case ComplexInt:    memcpy(val_ptr, &cival, size); break;
		case ComplexFloat:  memcpy(val_ptr, &cfval, size); break;
		case ComplexDouble: memcpy(val_ptr, &cdval, size); break;
		case Polar:         memcpy(val_ptr, &pval, size);  break;
		default: ;
	}
	
	return(val_ptr);
}

/************************************************************************
@Function: img_to_datatype
@Description:
	Convert an image to any new data type.
@Algorithm:
	This function is used in optional command-line arguments to indicate 
	a new data type for an image.
@Arguments:
	char letter			letter indicating data type.
	DataType newtype	data type (defined in rwimg.h).
@Returns:
	int					error code.
**************************************************************************/
int			img_to_datatype(Bimage* p, DataType newtype)
{
	if ( p->datatype == newtype ) return(0);
	
	int			err = 0;
	
	switch ( newtype ) {
		case UChar: 		err = img_to_byte(p);  break;
		case SChar: 		err = img_to_signed_byte(p);  break;
		case UShort: 		err = img_to_unsigned_short(p); break;
		case Short: 		err = img_to_short(p); break;
		case Int: 			err = img_to_int(p); 	break;
		case Float: 		err = img_to_float(p); break;
		case ComplexShort:  err = img_to_complex_short(p); break;
		case ComplexInt: 	err = img_to_complex_int(p); 	break;
		case ComplexFloat:  err = img_to_complex_float(p); break;
		case Polar: 		err = img_to_polar(p); break;
		default:			err = -1; break;
	}
	
	return(err);
}

/************************************************************************
@Function: img_to_bit
@Description:
	Converts an image to a bitmap.
@Algorithm:
	The image is thresholded and the bits written into a bitmap.
	The threshold is taken from the minimum and maximum:
		threshold = (min + max)/2
	Complex data types are converted to polar form and the intensities thresholded.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_to_bit(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->datatype >= ComplexShort ) img_complex2intensities(p);

    if ( p->colormodel == Bit ) return(-1);
	
    if ( p->colormodel == RGB ) img_RGB2gray(p);
	
    if ( p->colormodel == Index ) img_indexed2gray(p);
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    
    unsigned long		i, x, nyz = (unsigned long) p->y*p->z*p->c*p->n;
    unsigned long		datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
	unsigned long		new_x = (p->x-1)/8 + 1;
	unsigned long		bn = new_x*p->y*p->z*p->n;
    
    unsigned char* 		bdata = (unsigned char *) balloc(bn*sizeof(unsigned char));
    float				threshold = (p->max + p->min)/2;
	
    switch ( p->datatype ) {
    	case UChar:
			if ( verbose & VERB_LABEL )
    		    printf("Converting: unsigned char -> bitmap\n");
			for ( i=0; i<nyz; i++ )
				for ( x=0; x<p->x; x++ )
					if ( udata[i*p->x+x] > threshold ) bdata[i*new_x+x/8] |= 0x80 >> x%8;
			break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> bitmap\n");
			for ( i=0; i<nyz; i++ )
				for ( x=0; x<p->x; x++ )
					if ( cdata[i*p->x+x] > threshold ) bdata[i*new_x+x/8] |= 0x80 >> x%8;
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> bitmap\n");
			for ( i=0; i<nyz; i++ )
				for ( x=0; x<p->x; x++ )
					if ( sdata[i*p->x+x] > threshold ) bdata[i*new_x+x/8] |= 0x80 >> x%8;
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: unsigned short -> bitmap\n");
			for ( i=0; i<nyz; i++ )
				for ( x=0; x<p->x; x++ )
					if ( usdata[i*p->x+x] > threshold ) bdata[i*new_x+x/8] |= 0x80 >> x%8;
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: int -> bitmap\n");
			for ( i=0; i<nyz; i++ )
				for ( x=0; x<p->x; x++ )
					if ( idata[i*p->x+x] > threshold ) bdata[i*new_x+x/8] |= 0x80 >> x%8;
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> bitmap\n");
			for ( i=0; i<nyz; i++ )
				for ( x=0; x<p->x; x++ )
					if ( fdata[i*p->x+x] > threshold ) bdata[i*new_x+x/8] |= 0x80 >> x%8;
    	    break;
    	default: break;
    }
    
	if ( verbose & VERB_PROCESS )
    	printf("Threshold:                      %g\n\n", threshold);

	p->x = new_x*8;
	p->data = (char *) bdata;
	p->colormodel = Bit;
    p->datatype = UChar;
	bfree(udata, datasize*elementsize);
		
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_bit_to_byte
@Description:
	Converts an image from a bitmap to a byte type.
@Algorithm:
	The image is rescaled to fit into the range for a single byte (0 - 255).
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_bit_to_byte(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->datatype != UChar ) return(-1);

    if ( p->colormodel != Bit ) return(-1);
	    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
    
    unsigned char* 		bdata = (unsigned char *) p->data;	// Bitmap
    unsigned char* 		udata = (unsigned char *) balloc(datasize*sizeof(unsigned char));
	
	if ( verbose & VERB_LABEL )
    	printf("Converting: bitmap -> unsigned char\n");
	
	for ( i=0; i<datasize; i++ )
		if ( 0x80 & ( bdata[i/8] << i%8 ) )
			udata[i] = 255;

	p->data = (char *) udata;
	bfree(bdata, datasize/8);
	p->colormodel = Gray;
    p->datatype = UChar;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_byte
@Description:
	Converts an image to character or byte type.
@Algorithm:
	The image is rescaled to fit into the range for a single byte,
	with truncation of the data below 0 and above 255:
		new_data = data*255/(max-min)
	The minimum and maximum in the image structure determines how conversion
	is done.
	Complex data types are converted to polar form and the intensities returned.
	Signed byte images are shift up by 128 only if the minimum is less than zero.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_to_byte(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->datatype >= ComplexShort ) {
		img_complex2intensities(p);
		img_stats(p);
	}

    if ( p->datatype == UChar ) {
		p->scale = 1;
		return(-1);
	}
	
	if ( fabs(p->max - p->min) < 1e-37 ) p->max = p->min + 1;
	
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
	double				value;
    
    unsigned char* 		udata = (unsigned char *) balloc(datasize*sizeof(unsigned char));
    p->scale = 255/(p->max - p->min);
	
    switch ( p->datatype ) {
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> unsigned char\n");
    	    if ( p->min < 0 )
    	    	 for ( i=0; i<datasize; i++ ) udata[i] = cdata[i] + 128;
			else
				 for ( i=0; i<datasize; i++ ) udata[i] = cdata[i];
			p->scale = 1;
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> unsigned char\n");
    	    for ( i=0; i<datasize; i++ ) {
    	    	value = floor(p->scale*(sdata[i]-p->min)+0.5);
				if ( value > 255 ) udata[i] = 255;
				else if ( value < 0 ) udata[i] = 0;
				else udata[i] = (unsigned char) value;
			}
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: unsigned short -> unsigned char\n");
    	    for ( i=0; i<datasize; i++ ) {
    	    	value = floor(p->scale*(usdata[i]-p->min)+0.5);
				if ( value > 255 ) udata[i] = 255;
				else if ( value < 0 ) udata[i] = 0;
				else udata[i] = (unsigned char) value;
			}
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: int -> unsigned char\n");
    	    for ( i=0; i<datasize; i++ ) {
    	    	value = floor(p->scale*(idata[i]-p->min)+0.5);
				if ( value > 255 ) udata[i] = 255;
				else if ( value < 0 ) udata[i] = 0;
				else udata[i] = (unsigned char) value;
			}
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> unsigned char\n");
    	    for ( i=0; i<datasize; i++ ) {
    	    	value = floor(p->scale*(fdata[i]-p->min)+0.5);
				if ( value > 255 ) udata[i] = 255;
				else if ( value < 0 ) udata[i] = 0;
				else udata[i] = (unsigned char) value;
			}
    	    break;
    	default: break;
    }
    
	if ( verbose & VERB_PROCESS )
    	printf("Shift and scale:                %g %g\n\n", -p->min, p->scale);

	p->data = (char *) udata;
	bfree(cdata, datasize*elementsize);
	
    p->datatype = UChar;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_signed_byte
@Description:
	Converts a byte image to signed char.
@Algorithm:
	The input image is first converted to a byte image.
	If the range is witin 0 - 127, no shift is applied.
	Otherwise, the data is shifted by -128.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_to_signed_byte(Bimage* p)
{
	img_to_byte(p);
	
	int					shift = 128;
	
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;

	if ( verbose & VERB_LABEL )
		printf("Converting: unsigned char -> signed char\n");

	if ( p->max > 127 ) {
		for ( i=0; i<datasize; i++ )
			cdata[i] = udata[i] - shift;
		p->avg -= shift;
		p->min -= shift;
		p->max -= shift;
	}

    p->datatype = SChar;
	
	return(0);
}

/************************************************************************
@Function: img_to_unsigned_short
@Description:
	Converts an image to signed short or two-byte type.
@Algorithm:
	The image is rescaled to fit into the range for a signed short,
	(i.e., -32768 to 32767):
		new_data = data*scale/(max-min)
	The scale is 1 when the data is in the range, -32000 to 32000, or
	otherwise:
		scale = 32000/(max-min)
	The minimum and maximum in the image structure determines how conversion
	is done. If min and max are less than an integer apart, then:
		max = min + 1
	Complex data types are converted to polar form and the intensities returned.
	Colour images are converted to gray scale.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_to_unsigned_short(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
    if ( p->datatype == UShort ) return(-1);
    
	if ( p->datatype >= ComplexShort ) {
		error_show("Error: Complex unsigned short data type not supported", __FILE__, __LINE__);
		return(-1);
	}

	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    if ( p->colormodel == Index ) img_indexed2gray(p);
	
	if ( fabs(p->max - p->min) < 1 ) p->max = p->min + 1;
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    short*				sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
    
    unsigned short* 	usdata = (unsigned short *) balloc(datasize*sizeof(unsigned short));
	
	float				offset = 0;
	
	p->scale = 1;
		
	// Shrink the range if it is too large for an integer
	if ( p->min < 0 && p->min >= SHRT_MIN && p->max <= SHRT_MAX ) {
		offset = SHRT_MAX;
	} else {
    	p->scale = USHRT_MAX/(p->max - p->min);
		offset = -p->min;
	}
	
    switch ( p->datatype ) {
    	case UChar:
			if ( p->colormodel == Bit ) {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: bitmap -> short\n");
				for ( i=0; i<datasize; i++ )
					if ( 0x80 & ( udata[i/8] << i%8 ) )
						usdata[i] = 255;
			} else {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: unsigned char -> unsigned short\n");
    		    for ( i=0; i<datasize; i++ ) usdata[i] = udata[i];
			}
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> unsigned short\n");
			if ( p->min < 0 && p->max < 128 ) offset = 128;
			p->scale = 1;
    	    for ( i=0; i<datasize; i++ ) usdata[i] = (unsigned short)(cdata[i] + offset);
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: short -> unsigned short\n");
			offset = 0;
			if ( p->max > SHRT_MAX ) offset = SHRT_MIN;
			p->scale = 1;
    	    for ( i=0; i<datasize; i++ )
    	    	usdata[i] = (unsigned short) floor(sdata[i] + offset);
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: int -> unsigned short\n");
    	    for ( i=0; i<datasize; i++ )
    	    	usdata[i] = (unsigned short) floor(p->scale*(idata[i] + offset)+0.5);
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> unsigned short\n");
    	    for ( i=0; i<datasize; i++ )
    	    	usdata[i] = (unsigned short) floor(p->scale*(fdata[i] + offset)+0.5);
    	    break;
    	default: break;
    }
    
	if ( verbose & VERB_PROCESS )
    	printf("Shift and scale:                %g %g\n\n", offset, p->scale);

	p->data = (char *) usdata;
	if ( p->colormodel == Bit ) {
		bfree(udata, (datasize*elementsize)/8);
		p->colormodel = Gray;
	} else
		bfree(udata, datasize*elementsize);
	
	p->datatype = UShort;
		
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_short
@Description:
	Converts an image to signed short or two-byte type.
@Algorithm:
	The image is rescaled to fit into the range for a signed short,
	(i.e., -32768 to 32767):
		new_data = data*scale/(max-min)
	The scale is 1 when the data is in the range, -32000 to 32000, or
	otherwise:
		scale = 32000/(max-min)
	The minimum and maximum in the image structure determines how conversion
	is done. If min and max are less than an integer apart, then:
		max = min + 1
	Complex data types are converted to polar form and the intensities returned.
	Colour images are converted to gray scale.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_to_short(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
    if ( p->datatype == Short ) return(-1);
    
	if ( p->datatype >= ComplexShort ) img_complex2intensities(p);

	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    if ( p->colormodel == Index ) img_indexed2gray(p);
	
	if ( fabs(p->max - p->min) < 1 ) p->max = p->min + 1;
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    int* 	    		idata = (int *) p->data;
    float*  	    	fdata = (float *) p->data;
    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
    
    short* 	    		sdata = (short *) balloc(datasize*sizeof(short));
	
	float				offset = 0;
	
	p->scale = 1;
		
	// Shrink the range if it is too large for an integer
	if ( fabs(p->max) > 32767 || fabs(p->min) > 32767 ) {
    	p->scale = 32767/(p->max - p->min);
		offset = -(p->max + p->min)/2;
	}
	
	if ( p->max - p->min < 100 )
    	p->scale = 100/(p->max - p->min);
	
    switch ( p->datatype ) {
    	case UChar:
			if ( p->colormodel == Bit ) {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: bitmap -> short\n");
				for ( i=0; i<datasize; i++ )
					if ( 0x80 & ( udata[i/8] << i%8 ) )
						sdata[i] = 255;
			} else {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: unsigned char -> short\n");
    		    for ( i=0; i<datasize; i++ ) sdata[i] = udata[i];
			}
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> short\n");
    	    for ( i=0; i<datasize; i++ ) sdata[i] = cdata[i];
    	    break;
    	case Short:
    	case ComplexShort:
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: unsigned short -> short\n");
			if ( p->max > 32767 ) offset = -32768;
			p->scale = 1;
    	    for ( i=0; i<datasize; i++ )
    	    	sdata[i] = (short) floor(1.0*usdata[i] + offset);
    	    break;
    	case Int:
    	case ComplexInt:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: int -> short\n");
    	    for ( i=0; i<datasize; i++ )
    	    	sdata[i] = (short) floor(p->scale*(idata[i] + offset)+0.5);
    	    break;
    	case Float:
    	case ComplexFloat:
    	case Polar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> short\n");
    	    for ( i=0; i<datasize; i++ )
    	    	sdata[i] = (short) floor(p->scale*(fdata[i] + offset)+0.5);
    	    break;
    	default: break;
    }
    
	if ( verbose & VERB_PROCESS )
    	printf("Shift and scale:                %g %g\n\n", offset, p->scale);

	p->data = (char *) sdata;
	if ( p->colormodel == Bit ) {
		bfree(udata, (datasize*elementsize)/8);
		p->colormodel = Gray;
	} else
		bfree(udata, datasize*elementsize);
	
	if ( p->datatype < ComplexShort )
    	p->datatype = Short;
	else
    	p->datatype = ComplexShort;
		
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_int
@Description:
	Converts an image to integer or four-byte type.
@Algorithm:
	The image is rescaled to fit into the range for a signed integer
	(i.e., -256*256*256*128 to 256*256*256*128-1):
		new_data = data*scale/(max-min)
	The scale is 1 when the data is in the range, -2e9 to 2e9, or
	otherwise:
		scale = 2e9/(max-min)
	The minimum and maximum in the image structure determines how conversion
	is done. If min and max are less than an integer apart, then:
		max = min + 1
	Complex data types are converted to polar form and the intensities returned.
	Colour images are converted to gray scale.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_to_int(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->datatype >= ComplexShort ) img_complex2intensities(p);

	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    if ( p->colormodel == Index ) img_indexed2gray(p);
	
    if ( p->datatype == Int ) return(-1);
    
	if ( fabs(p->max - p->min) < 1 ) p->max = p->min + 1;
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    float*  	    	fdata = (float *) p->data;
    
    unsigned long		i, datasize = (unsigned long) p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
    
    int* 	    		idata = (int *) balloc(datasize*sizeof(int));
	
	p->scale = 1;
		
	// Shrink the range if it is too large for an integer
	if ( p->max > 2e9 || p->min < -2e9 )
    	p->scale = 2e9/(p->max - p->min);
	
    switch ( p->datatype ) {
    	case UChar:
			if ( p->colormodel == Bit ) {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: bitmap -> int\n");
				for ( i=0; i<datasize; i++ )
					if ( 0x80 & ( udata[i/8] << i%8 ) )
						idata[i] = 255;
			} else {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: unsigned char -> int\n");
    		    for ( i=0; i<datasize; i++ ) idata[i] = udata[i];
			}
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> int\n");
    	    for ( i=0; i<datasize; i++ ) idata[i] = cdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> int\n");
    	    for ( i=0; i<datasize; i++ ) idata[i] = sdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
	    	    printf("Converting: unsigned short -> int\n");
    	    for ( i=0; i<datasize; i++ )
    	    	idata[i] = (int) floor(p->scale*usdata[i]+0.5);
    	    break;
    	case Float:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: float -> int\n");
    	    for ( i=0; i<datasize; i++ )
    	    	idata[i] = (int) floor(p->scale*fdata[i]+0.5);
    	    break;
    	default: break;
    }
    
	if ( verbose & VERB_PROCESS )
    	printf("Scale:                          %g\n\n", p->scale);

	p->data = (char *) idata;
	if ( p->colormodel == Bit ) {
		bfree(udata, (datasize*elementsize)/8);
		p->colormodel = Gray;
	} else
		bfree(udata, datasize*elementsize);
	
	if ( p->datatype < ComplexShort )
    	p->datatype = Int;
	else
    	p->datatype = ComplexInt;

	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_to_float
@Description:
	Converts an image to four-byte floating point type.
@Algorithm:
	The data is converted to floating point without rescaling.
	Complex data types are converted to polar form and the intensities returned.
	Colour images are converted to gray scale.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_to_float(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->datatype >= ComplexShort ) img_complex2intensities(p);

    if ( p->datatype == Float ) return(-1);
    
	if ( p->colormodel == RGB ) img_RGB2gray(p);
    
    if ( p->colormodel == Index ) img_indexed2gray(p);
	
    unsigned char* 		udata = (unsigned char *) p->data;
    signed char* 		cdata = (signed char *) p->data;
    unsigned short* 	usdata = (unsigned short *) p->data;
    short* 	    		sdata = (short *) p->data;
    int* 	    		idata = (int *) p->data;
    
    unsigned long		i, datasize = p->x*p->y*p->z*p->c*p->n;
	unsigned long		elementsize = gettypesize(p->datatype);
    float*  	    	fdata = (float *) balloc(datasize*sizeof(float));
    
    switch ( p->datatype ) {
    	case UChar:
			if ( p->colormodel == Bit ) {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: bitmap -> float\n");
				for ( i=0; i<datasize; i++ )
					if ( 0x80 & ( udata[i/8] << i%8 ) )
						fdata[i] = 1;
			} else {
				if ( verbose & VERB_LABEL )
    		    	printf("Converting: unsigned char -> float\n\n");
    		    for ( i=0; i<datasize; i++ ) fdata[i] = udata[i];
			}
    	    break;
    	case SChar:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: signed char -> float\n\n");
    	    for ( i=0; i<datasize; i++ ) fdata[i] = cdata[i];
    	    break;
    	case UShort:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: unsigned short -> float\n\n");
    	    for ( i=0; i<datasize; i++ ) fdata[i] = usdata[i];
    	    break;
    	case Short:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: short -> float\n\n");
    	    for ( i=0; i<datasize; i++ ) fdata[i] = sdata[i];
    	    break;
    	case Int:
			if ( verbose & VERB_LABEL )
    	    	printf("Converting: int -> float\n\n");
    	    for ( i=0; i<datasize; i++ ) fdata[i] = idata[i];
    	    break;
    	default: break;
    }
    
	p->data = (char *) fdata;
	if ( p->colormodel == Bit ) {
		bfree(udata, (datasize*elementsize)/8);
		p->colormodel = Gray;
	} else
		bfree(udata, datasize*elementsize);
	
	if ( p->datatype < ComplexShort )
    	p->datatype = Float;
	else
    	p->datatype = ComplexFloat;
	
	img_stats(p);
	
	return(0);
}

/************************************************************************
@Function: img_gray2RGB
@Description:
	Converts an image from gray scale to RGB (red-green-blue).
@Algorithm:
	The image is first converted to byte format.
	Each byte is then expanded to three with the same value.
	The key variable that is changed is the number channels, going from
	one for gray scale images to three for RGB.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int 		img_gray2RGB(Bimage* p)
{
	return(img_gray2color(p, 0));
}

/************************************************************************
@Function: img_gray2RGBA
@Description:
	Converts an image from gray scale to RGB (red-green-blue-alpha).
@Algorithm:
	The image is first converted to byte format.
	Each byte is then expanded to three with the same value.
	The key variable that is changed is the number channels, going from
	one for gray scale images to three for RGB.
	The fourth channel (transparency or alpha) is set to 0 if the
	input is zero and 255 otherwise.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int 		img_gray2RGBA(Bimage* p)
{
	return(img_gray2color(p, 1));
}

/************************************************************************
@Function: img_gray2color
@Description:
	Converts an image from gray scale to RGB or RGBA.
@Algorithm:
	The image is first converted to byte format.
	Each byte is then expanded to three with the same value.
	The key variable that is changed is the number channels, going from
	one for gray scale images to three for RGB.
	If the alpha flag is nonzero the fourth channel (transparency or alpha) 
	is set to 0 if the input is zero and 255 otherwise.
@Arguments:
	Bimage* p		image (replaced by new image).
	int alpha		flag to indicate alpha channel.
@Returns:
	int				error code.
**************************************************************************/
int 		img_gray2color(Bimage* p, int alpha)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->colormodel >= RGB ) return(0);
	
	img_to_byte(p);
	
	unsigned long   i, j, c, nc = 3, datasize = (unsigned long) p->x*p->y*p->z*p->n;
	unsigned char   max = 255;
	if ( p->max < 255 ) max = (unsigned char) p->max;
	if ( alpha ) nc = 4;
	
	char*		data = (char *) p->data;
	char*		rgbdata = (char *) balloc(nc*datasize*sizeof(char));
	
	for ( i=0; i<datasize; i++ ) {
		for ( c=0, j=nc*i; c<3; c++, j++ )
			rgbdata[j] = data[i];
		if ( alpha )
			if ( data[i] ) rgbdata[j] = max;
	}
	
	if ( verbose & VERB_LABEL )  {
		if ( alpha ) printf("Converting: unsigned char -> RGBA:\n");
		else printf("Converting: unsigned char -> RGB:\n");
		printf("Scale:                          %g\n\n", p->scale);
	}
		
	p->data = (char *) rgbdata;
	bfree(data, datasize*sizeof(char));
	
	p->datatype = UChar;
	p->colormodel = RGB;
	if ( alpha ) p->colormodel = RGBA;
	p->c = nc;
	
	return(0);
}

/************************************************************************
@Function: img_RGB2gray
@Description:
	Converts an image from RGB (red-green-blue) to gray scale.
@Algorithm:
	The key variable that is changed is the number channels, going from
	three for RGB images to one for gray scale.
	The new gray scale value is just the average of the three RGB values.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int 		img_RGB2gray(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->colormodel != RGB && p->colormodel != RGBA ) return(0);
	
	if ( p->datatype > Int ) img_to_byte(p);
	
	if ( verbose & VERB_LABEL )  {
		printf("Converting: RGB -> gray scale\n");
		printf("\n");
	}
		
	unsigned long   i, j, c, nc = 3, datasize = (unsigned long) p->x*p->y*p->z*p->n;
	unsigned long   datatypesize = gettypesize(p->datatype);

	unsigned char*	ucrgb = (unsigned char *) p->data;
	unsigned short*	usrgb = (unsigned short *) p->data;
	int*			irgb = (int *) p->data;
	
	unsigned char*	ucdata = (unsigned char *) balloc(datasize*datatypesize);
	unsigned short*	usdata = (unsigned short *) ucdata;
	int*			idata = (int *) ucdata;
	
	switch ( p->datatype ) {
		case UChar:
		case SChar:
			for ( i=0; i<datasize; i++ )
				for ( c=0, j=p->c*i; c<nc; c++, j++ )
					ucdata[i] += ucrgb[j]/nc;
			break;
		case UShort:
		case Short:
			for ( i=0; i<datasize; i++ )
				for ( c=0, j=p->c*i; c<nc; c++, j++ )
					usdata[i] = (short) (usrgb[j]/nc);
			break;
		case Int:
			for ( i=0; i<datasize; i++ )
				for ( c=0, j=p->c*i; c<nc; c++, j++ )
					idata[i] = (int) (irgb[j]/nc);
			break;
	}
	
	p->data = (char *) ucdata;
	bfree(ucrgb, datasize*p->c*datatypesize);
	
	p->colormodel = Gray;
	p->c = 1;
	p->scale = 1;
	
	return(0);
}

/************************************************************************
@Function: img_indexed2gray
@Description:
	Converts an image from indexed colour map to RGB (red-green-blue).
@Algorithm:
	The RGB values of the colour map is applied to each pixel.
	The key variable that is changed is the number channels, going from
	one for gray scale images to three for RGB.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int 		img_indexed2gray(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->colormodel != Index || p->colormap == NULL ) return(0);
	
	unsigned long   i, j, datasize = (unsigned long) p->x*p->y*p->z*p->n;
	unsigned char*	data = (unsigned char *) p->data;
	
	for ( i=0; i<datasize; i++ ) {
		j = 3*data[i];
		data[i] = (p->colormap[j] + p->colormap[j+1] + p->colormap[j+2])/3;
	}
	
	if ( verbose & VERB_LABEL )
		printf("Converting: indexed -> gray scale\n");
		
	p->datatype = UChar;
	p->colormodel = Gray;
	
	return(0);
}

/************************************************************************
@Function: img_indexed2RGB
@Description:
	Converts an image from indexed colour map to RGB (red-green-blue).
@Algorithm:
	The RGB values of the colour map is applied to each pixel.
	The key variable that is changed is the number channels, going from
	one for gray scale images to three for RGB.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int 		img_indexed2RGB(Bimage* p)
{
	if ( p->dataflag < 1 ) return(-1);
	
	if ( p->colormodel != Index || p->colormap == NULL ) return(0);
	
	unsigned long   i, j, datasize = (unsigned long) p->x*p->y*p->z*p->n;
	char*			data = (char *) p->data;
	char*			rgbdata = (char *) balloc(3*datasize*sizeof(char));
	
	for ( i=0; i<datasize; i++ ) {
		j = 3*data[i];
		rgbdata[3*i] = p->colormap[j];
		rgbdata[3*i+1] = p->colormap[j+1];
		rgbdata[3*i+2] = p->colormap[j+2];
	}
	
	if ( verbose & VERB_LABEL )
		printf("Converting: indexed -> RGB\n");
		
	p->data = (char *) rgbdata;
	bfree(data, datasize*sizeof(char));
	
	p->datatype = UChar;
	p->colormodel = RGB;
	p->c = 3;
	
	return(0);
}

/************************************************************************
@Function: img_add_alpha
@Description:
	Adds an alpha channel to generate a RGBA (red-green-blue-alpha) image.
@Algorithm:
	The new alpha channel contains a zero for original black and 255
	for original nonblack RGB colours.
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_add_alpha(Bimage* p)
{
	if ( p->colormodel < RGB ) return(img_gray2RGBA(p));
	
	if ( p->colormodel != RGB ) return(0);
	
	unsigned long   i, j, k, c, s, datasize = (unsigned long) p->x*p->y*p->z*p->n;
	unsigned char*  rgb = (unsigned char *) p->data;
	unsigned char*  rgba = (unsigned char *) balloc(4*datasize*sizeof(unsigned char));
	
	if ( verbose & VERB_LABEL )
		printf("Converting: RGB -> RGBA\n");
		
	for ( i=0; i<datasize; i++ ) {
		for ( c=0, j=3*i, k=4*i, s=0; c<3; c++, j++, k++ ) {
			rgba[k] = rgb[j];
			s += rgb[j];
		}
		if ( s ) rgba[k] = 255;
	}
	
	bfree(rgb, 3*datasize*sizeof(char));
	
	p->data = (char *) rgba;
	p->c = 4;
	p->colormodel = RGBA;
	
	return(0);
}

/************************************************************************
@Function: img_remove_alpha
@Description:
	Removes an alpha channel from a RGBA (red-green-blue-alpha) image.
@Algorithm:
	The new image is just RGB (red-green-blue).
@Arguments:
	Bimage* p		image (replaced by new image).
@Returns:
	int				error code.
**************************************************************************/
int			img_remove_alpha(Bimage* p)
{
	if ( p->colormodel != RGBA ) return(0);
	
	unsigned long   i, j, k, c, datasize = (unsigned long) p->x*p->y*p->z*p->n;
	char*			rgba = (char *) p->data;
	char*			rgb = (char *) balloc(3*datasize*sizeof(char));
	
	if ( verbose & VERB_LABEL )
		printf("Converting: RGBA -> RGB\n");
		
	for ( i=0; i<datasize; i++ ) {
		for ( c=0, j=3*i, k=4*i; c<3; c++, j++, k++ )
			rgb[j] = rgba[k];
	}
	
	bfree(rgba, 4*datasize*sizeof(char));
	
	p->data = (char *) rgb;
	p->c = 3;
	p->colormodel = RGB;
	
	return(0);
}


