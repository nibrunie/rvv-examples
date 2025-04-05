#ifndef _AES_GCM_VECTORS
#define _AES_GCM_VECTORS

#include "gcmEncryptExtIV128.h"
#include "gcmDecrypt128.h"
#include "gcmDecrypt256.h"
#include "gcmEncryptExtIV256.h"

static const struct aes_gcm_test_suite gcm_suites[] = {
    {.name = "gcmEncryptExtIV128", .count = 7875, .keylen = 128, .tests = gcmEncryptExtIV128},
    {.name = "gcmDecrypt128", .count = 7875, .keylen = 128, .tests = gcmDecrypt128},
    {.name = "gcmDecrypt256", .count = 7875, .keylen = 256, .tests = gcmDecrypt256},
    {.name = "gcmEncryptExtIV256", .count = 7875, .keylen = 256, .tests = gcmEncryptExtIV256},
};
#endif
