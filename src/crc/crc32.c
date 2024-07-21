#include <stdint.h>
#include <stdlib.h>

typedef uint32_t u32;

// number of message bits processed during each iteration of the CRC
#ifndef CRC_LE_BITS
#  define CRC_LE_BITS 1
#endif

#ifndef CRC_BE_BITS
#  define CRC_BE_BITS 1
#endif

/** copied from
 * https://elixir.bootlin.com/linux/v4.4/source/lib/gen_crc32table.c#L59
*/

#define ENTRIES_PER_LINE 4

#if CRC_LE_BITS > 8
# define LE_TABLE_ROWS (CRC_LE_BITS/8)
# define LE_TABLE_SIZE 256
#else
# define LE_TABLE_ROWS 1
# define LE_TABLE_SIZE (1 << CRC_LE_BITS)
#endif

#if CRC_BE_BITS > 8
# define BE_TABLE_ROWS (CRC_BE_BITS/8)
# define BE_TABLE_SIZE 256
#else
# define BE_TABLE_ROWS 1
# define BE_TABLE_SIZE (1 << CRC_BE_BITS)
#endif

uint32_t crc32table_le[LE_TABLE_ROWS][256] = {0};
uint32_t crc32table_be[BE_TABLE_ROWS][256] = {0};

/**
 * crc32init_le() - allocate and initialize LE table data
 *
 * crc is the crc of the byte i; other entries are filled in based on the
 * fact that crctable[i^j] = crctable[i] ^ crctable[j].
 *
 */
void crc32init_le_generic(const uint32_t polynomial,
				 uint32_t (*tab)[256])
{
	unsigned i, j;
	uint32_t crc = 1;

	tab[0][0] = 0;

	for (i = LE_TABLE_SIZE >> 1; i; i >>= 1) {
		crc = (crc >> 1) ^ ((crc & 1) ? polynomial : 0);
		for (j = 0; j < LE_TABLE_SIZE; j += 2 * i)
			tab[0][i + j] = crc ^ tab[0][j];
	}
	for (i = 0; i < LE_TABLE_SIZE; i++) {
		crc = tab[0][i];
		for (j = 1; j < LE_TABLE_ROWS; j++) {
			crc = tab[0][crc & 0xff] ^ (crc >> 8);
			tab[j][i] = crc;
		}
	}
}

void crc32init_le(const uint32_t polynomial) {
    crc32init_le_generic(polynomial, crc32table_le);
}

/**
 * Copied from https://elixir.bootlin.com/linux/v4.4/source/lib/crc32.c#L145
 */


/**
 * crc32_le_generic() - Calculate bitwise little-endian Ethernet AUTODIN II
 *			CRC32/CRC32C
 * @crc: seed value for computation.  ~0 for Ethernet, sometimes 0 for other
 *	 uses, or the previous crc32/crc32c value if computing incrementally.
 * @p: pointer to buffer over which CRC32/CRC32C is run
 * @len: length of buffer @p
 * @polynomial: CRC32/CRC32c LE polynomial
 */
u32 crc32_le_generic(u32 crc, unsigned char const *p,
					 size_t len,
					 u32 polynomial)
{
    u32 (*tab)[256] = crc32table_le;
	int i;
    # if CRC_LE_BITS == 1
	while (len--) {
		crc ^= *p++;
		for (i = 0; i < 8; i++)
			crc = (crc >> 1) ^ ((crc & 1) ? polynomial : 0);
	}
    # elif CRC_LE_BITS == 2
	while (len--) {
		crc ^= *p++;
		crc = (crc >> 2) ^ tab[0][crc & 3];
		crc = (crc >> 2) ^ tab[0][crc & 3];
		crc = (crc >> 2) ^ tab[0][crc & 3];
		crc = (crc >> 2) ^ tab[0][crc & 3];
	}
# elif CRC_LE_BITS == 4
	while (len--) {
		crc ^= *p++;
		crc = (crc >> 4) ^ tab[0][crc & 15];
		crc = (crc >> 4) ^ tab[0][crc & 15];
	}
# elif CRC_LE_BITS == 8
	/* aka Sarwate algorithm */
	while (len--) {
		crc ^= *p++;
		crc = (crc >> 8) ^ tab[0][crc & 255];
	}
# else
#   error "unsupported CRC_LE_BITS value " ## CRC_LE_BITS
# endif
	return crc;
}


/**
 * crc32_be_generic() - Calculate bitwise big-endian Ethernet AUTODIN II CRC32
 * @crc: seed value for computation.  ~0 for Ethernet, sometimes 0 for
 *	other uses, or the previous crc32 value if computing incrementally.
 * @p: pointer to buffer over which CRC32 is run
 * @len: length of buffer @p
 * @tab: big-endian Ethernet table
 * @polynomial: CRC32 BE polynomial
 */
u32 crc32_be_generic(u32 crc, unsigned char const *p,
					  size_t len,
					  u32 polynomial)
{
u32 (*tab)[256] = crc32table_be;
#if CRC_BE_BITS == 1
	int i;
	while (len--) {
		crc ^= *p++ << 24;
		for (i = 0; i < 8; i++)
			crc =
			    (crc << 1) ^ ((crc & 0x80000000) ? polynomial :
					  0);
	}
# elif CRC_BE_BITS == 2
	while (len--) {
		crc ^= *p++ << 24;
		crc = (crc << 2) ^ tab[0][crc >> 30];
		crc = (crc << 2) ^ tab[0][crc >> 30];
		crc = (crc << 2) ^ tab[0][crc >> 30];
		crc = (crc << 2) ^ tab[0][crc >> 30];
	}
# elif CRC_BE_BITS == 4
	while (len--) {
		crc ^= *p++ << 24;
		crc = (crc << 4) ^ tab[0][crc >> 28];
		crc = (crc << 4) ^ tab[0][crc >> 28];
	}
# elif CRC_BE_BITS == 8
	while (len--) {
		crc ^= *p++ << 24;
		crc = (crc << 8) ^ tab[0][crc >> 24];
	}
# else
# endif
	return crc;
}

void crc32init_be_generic(const uint32_t polynomial,
				 uint32_t (*tab)[256])
{
		unsigned i, j;
	uint32_t crc = 0x80000000;

	crc32table_be[0][0] = 0;

	for (i = 1; i < BE_TABLE_SIZE; i <<= 1) {
		crc = (crc << 1) ^ ((crc & 0x80000000) ? polynomial : 0);
		for (j = 0; j < i; j++)
			crc32table_be[0][i + j] = crc ^ crc32table_be[0][j];
	}
	for (i = 0; i < BE_TABLE_SIZE; i++) {
		crc = crc32table_be[0][i];
		for (j = 1; j < BE_TABLE_ROWS; j++) {
			crc = crc32table_be[0][(crc >> 24) & 0xff] ^ (crc << 8);
			crc32table_be[j][i] = crc;
		}
	}
}

void crc32init_be(const uint32_t polynomial) {
    crc32init_be_generic(polynomial, crc32table_be);
}