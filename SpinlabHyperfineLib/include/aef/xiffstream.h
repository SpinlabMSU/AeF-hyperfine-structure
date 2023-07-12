#pragma once
/*
 */
#ifndef _AEF_XIFFSTREAM_H
#define _AEF_XIFFSTREAM_H 1
#include "aef.h"
#include <filesystem>
#include <istream>
#include <vector>


/*
* XIFF is a 64-bit extension of RIFF (resource interchange file format).
* used here because I already wrote this before
 */
namespace xiff {

/*
 * A traditional Four Character Code (4 ascii characters == 32 bits)
 * Note that these generally aren't zero-terminated, so don't try to use them as strings!
 */
union fourcc {
  uint8_t cc[4];
  uint32_t ucode;
};

// common codes
namespace common_cc {
constexpr fourcc xiff = fourcc({'X', 'I', 'F', 'F'});
constexpr fourcc fmt_ = fourcc({'f', 'm', 't', ' '});
constexpr fourcc list = fourcc({'l', 'i', 's', 't'});
constexpr fourcc grup = fourcc({'g', 'r', 'p', ' '});
}; // namespace common_cc

struct chunk_hdr {
  // defines xiff chunk type
  fourcc type;
  // chunk format version
  uint16_t version;
  // flags are chunk-format
  uint16_t flags;
  // length in bytes
  uint64_t length;
  // followed by data
};

struct xiff_hdr {
  chunk_hdr file_hdr;
  // file type
  fourcc ftype;
  // extra specifier
  fourcc extra;
  // data begins here
};

class xiffstream {
  void *dataptr;
};

} // namespace xiff

#endif // _AEF_XIFFSTREAM_H