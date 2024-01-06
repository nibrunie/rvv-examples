/** return the value of selected perf counter
 * 
 * perf counter is selected through a macro:
 * - defining COUNT_INSTRET selects the instret counter
 *    The instret counter counts the number of retired (executed) instructions.
 * - defining COUNT_CYCLE selects cycle count
*/
static unsigned long read_perf_counter(void)
{
  unsigned long counter_value;
#if defined(COUNT_INSTRET)
  asm volatile ("rdinstret %0" : "=r" (counter_value));
#elif defined(COUNT_CYCLE)
  asm volatile ("rdcycle %0" : "=r" (counter_value));
#else
  // instret is also the default
  asm volatile ("rdinstret %0" : "=r" (counter_value));
#endif
  return counter_value;
}