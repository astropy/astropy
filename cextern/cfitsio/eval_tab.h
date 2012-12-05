typedef union {
    int    Node;        /* Index of Node */
    double dbl;         /* real value    */
    long   lng;         /* integer value */
    char   log;         /* logical value */
    char   str[MAX_STRLEN];    /* string value  */
} FFSTYPE;
#define	BOOLEAN	258
#define	LONG	259
#define	DOUBLE	260
#define	STRING	261
#define	BITSTR	262
#define	FUNCTION	263
#define	BFUNCTION	264
#define	IFUNCTION	265
#define	GTIFILTER	266
#define	REGFILTER	267
#define	COLUMN	268
#define	BCOLUMN	269
#define	SCOLUMN	270
#define	BITCOL	271
#define	ROWREF	272
#define	NULLREF	273
#define	SNULLREF	274
#define	OR	275
#define	AND	276
#define	EQ	277
#define	NE	278
#define	GT	279
#define	LT	280
#define	LTE	281
#define	GTE	282
#define	POWER	283
#define	NOT	284
#define	INTCAST	285
#define	FLTCAST	286
#define	UMINUS	287
#define	ACCUM	288
#define	DIFF	289


extern FFSTYPE fflval;
