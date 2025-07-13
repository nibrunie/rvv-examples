

#if defined(LINUX_PERF_COUNT)
// Code copied from: https://learn.arm.com/learning-paths/servers-and-cloud-computing/arm_pmu/perf_event_open/

#include <linux/perf_event.h> /* Definition of PERF_* constants */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/syscall.h> /* Definition of SYS_* constants */
#include <unistd.h>
#include <inttypes.h>


// Executes perf_event_open syscall and makes sure it is successful or exit
static long perf_event_open(struct perf_event_attr *hw_event, pid_t pid, int cpu, int group_fd, unsigned long flags){
  int fd;
  fd = syscall(SYS_perf_event_open, hw_event, pid, cpu, group_fd, flags);
  if (fd == -1) {
    fprintf(stderr, "Error creating event");
    exit(EXIT_FAILURE);
  }

  return fd;
}


int perf_count_init() {
  int fd;
  struct perf_event_attr  pe;

  // Configure the event to count
  memset(&pe, 0, sizeof(struct perf_event_attr));
  pe.type = PERF_TYPE_HARDWARE;
  pe.size = sizeof(struct perf_event_attr);
  pe.config = PERF_COUNT_HW_CPU_CYCLES; // PERF_COUNT_HW_INSTRUCTIONS;
  pe.disabled = 1;
  pe.exclude_kernel = 1;   // Do not measure instructions executed in the kernel
  pe.exclude_hv = 1;  // Do not measure instructions executed in a hypervisor

  // Create the event
  fd = perf_event_open(&pe, 0, -1, -1, 0);

  //Reset counters and start counting
  ioctl(fd, PERF_EVENT_IOC_RESET, 0);
  // Example code to count through

  return fd;
}
void perf_count_start(int fd) {
  ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
}

void perf_count_stop(int fd) {
  // Stop counting
  ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
}


unsigned long read_perf_counter(int fd) {
  uint64_t  val;
  // Read and print result
  read(fd, &val, sizeof(val));
  return val;
}

void perf_count_cleanup(int fd) {
  // Clean up file descriptor
  close(fd);
}
#else // defined(LINUX_PERF_COUNT)

static inline int perf_count_init(void) { return -1;};
static inline void perf_count_start(int fd)  {};
static inline void perf_count_stop(int fd)   {};
static inline int perf_count_cleanup(int fd) {};

/** return the value of selected perf counter
 * 
 * perf counter is selected through a macro:
 * - defining COUNT_INSTRET selects the instret counter
 *    The instret counter counts the number of retired (executed) instructions.
 * - defining COUNT_CYCLE selects cycle count
*/
static unsigned long read_perf_counter(int fd)
{
  unsigned long counter_value;
#if defined(COUNT_INSTRET)
#define PERF_METRIC "instruction"
  asm volatile ("rdinstret %0" : "=r" (counter_value));
#elif defined(COUNT_CYCLE)
#define PERF_METRIC "cycle"
  asm volatile ("rdcycle %0" : "=r" (counter_value));
#else
  // instret is also the default
#define PERF_METRIC "instruction"
  asm volatile ("rdinstret %0" : "=r" (counter_value));
#endif
  return counter_value;
}
#endif // defined(LINUX_PERF_COUNT)