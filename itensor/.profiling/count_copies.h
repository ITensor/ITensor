#ifndef __ITENSOR_COUNT_COPIES_H
#define __ITENSOR_COUNT_COPIES_H

//#define COUNT_COPIES

#ifdef COUNT_COPIES
#define IF_COUNT_COPIES(x) { x }
#else
#define IF_COUNT_COPIES(x) { }
#endif

#ifdef COUNT_COPIES

extern int copycount;
#ifdef THIS_IS_MAIN
int copycount = 0;
#endif //THIS_IS_MAIN

#endif

#endif
