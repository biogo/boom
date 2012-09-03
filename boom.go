// Copyright ©2012 Dan Kortschak <dan.kortschak@adelaide.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// Package boom is a wrapper for the samtools bam library.
package boom

// http://samtools.sourceforge.net/samtools/sam/index.html

/*
#cgo CFLAGS: -g -O2 -fPIC -m64 -pthread
#cgo LDFLAGS: -lz
#include "sam.h"
#include "bam_endian.h"
void bam_init_header_hash(bam_header_t *header);
void bam_destroy_header_hash(bam_header_t *header);
void setBin(bam1_t *b, uint16_t bin)        { b->core.bin = bin; }
void setQual(bam1_t *b, uint8_t flag)       { b->core.flag = flag; }
void setLQname(bam1_t *b, uint8_t l_qname)  { b->core.l_qname = l_qname; }
void setFlag(bam1_t *b, uint16_t flag)      { b->core.flag = flag; }
void setNCigar(bam1_t *b, uint16_t n_cigar) { b->core.n_cigar = n_cigar; }
*/
import "C"

import (
	"encoding/binary"
	"errors"
	"fmt"
	"io"
	"reflect"
	"runtime"
	"unsafe"
)

var (
	valueIsNil       = fmt.Errorf("boom: value is nil")
	notBamFile       = fmt.Errorf("boom: not bam file")
	couldNotAllocate = fmt.Errorf("boom: could not allocate")
	cannotAddr       = fmt.Errorf("boom: cannot address value")
	bamIsBigEndian   = C.bam_is_big_endian() == 1
	endian           = [2]binary.ByteOrder{
		binary.LittleEndian,
		binary.BigEndian,
	}[C.bam_is_big_endian()]
)

var (
	noHeader = errors.New("boom: no header")
)

// Verbosity sets and returns the level of debugging information emitted on stderr by libbam.
// The level of verbosity intrepreted by libbam ranges from 0 to 3 inclusive, with lower values
// being less verbose. Passing values of v outside this range do not alter verbosity.
func Verbosity(v int) int {
	if 0 <= v && v <= 3 {
		C.bam_verbose = C.int(v)
	}
	return int(C.bam_verbose)
}

// A bamRecord wraps the bam1_t BAM record.
type bamRecord struct {
	b *C.bam1_t
}

// newBamRecord creates a new bamRecord wrapping b or a newly malloc'd bam1_t if b is nil,
// and setting a finaliser that C.free()s the contained bam1_t.
// newBamRecord should always be used unless the bamRecord will be explicitly memory managed, or
// wraps a bam1_t that will be memory managed elsewhere.
func newBamRecord(b *C.bam1_t) (br *bamRecord, err error) {
	if b == nil {
		b = (*C.bam1_t)(unsafe.Pointer(C.malloc((C.size_t)(unsafe.Sizeof(C.bam1_t{})))))

		if b == nil {
			return nil, couldNotAllocate
		}
		*b = C.bam1_t{}
	}

	br = &bamRecord{b}
	runtime.SetFinalizer(br, (*bamRecord).bamRecordFree)

	return
}

// The following methods are helpers to safely return bam1_t field values.
// All first check that the pointer to the bam1_t is not nil and convert to the appropriate
// Go type.
func (br *bamRecord) tid() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.core.tid)
}
func (br *bamRecord) setTid(tid int32) {
	if br.b == nil {
		panic(valueIsNil)
	}
	br.b.core.tid = C.int32_t(tid)
}
func (br *bamRecord) pos() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.core.pos)
}
func (br *bamRecord) setPos(pos int32) {
	if br.b == nil {
		panic(valueIsNil)
	}
	br.b.core.pos = C.int32_t(pos)
}
func (br *bamRecord) bin() uint16 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return uint16(br.b.core.bin)
}
func (br *bamRecord) setBin(bin uint16) {
	if br.b == nil {
		panic(valueIsNil)
	}
	C.setBin(br.b, C.uint16_t(bin))
}
func (br *bamRecord) qual() byte {
	if br.b == nil {
		panic(valueIsNil)
	}
	return byte(br.b.core.qual)
}
func (br *bamRecord) setQual(qual byte) {
	if br.b == nil {
		panic(valueIsNil)
	}
	C.setQual(br.b, C.uint8_t(qual))
}
func (br *bamRecord) lQname() byte {
	if br.b == nil {
		panic(valueIsNil)
	}
	return byte(br.b.core.l_qname)
}
func (br *bamRecord) setLQname(lQname byte) {
	if br.b == nil {
		panic(valueIsNil)
	}
	C.setLQname(br.b, C.uint8_t(lQname))
}
func (br *bamRecord) flag() Flags {
	if br.b == nil {
		panic(valueIsNil)
	}
	return Flags(br.b.core.flag)
}
func (br *bamRecord) setFlag(flags Flags) {
	if br.b == nil {
		panic(valueIsNil)
	}
	C.setFlag(br.b, C.uint16_t(flags))
}
func (br *bamRecord) nCigar() uint16 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return uint16(br.b.core.n_cigar)
}
func (br *bamRecord) setNCigar(nCigar uint16) {
	if br.b == nil {
		panic(valueIsNil)
	}
	C.setNCigar(br.b, C.uint16_t(nCigar))
}
func (br *bamRecord) lQseq() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.core.l_qseq)
}
func (br *bamRecord) setLQseq(lQseq int32) {
	if br.b == nil {
		panic(valueIsNil)
	}
	br.b.core.l_qseq = C.int32_t(lQseq)
}
func (br *bamRecord) mtid() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.core.mtid)
}
func (br *bamRecord) setMtid() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.core.mtid)
}
func (br *bamRecord) mpos() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.core.mpos)
}
func (br *bamRecord) setMpos(mpos int32) {
	if br.b == nil {
		panic(valueIsNil)
	}
	br.b.core.mpos = C.int32_t(mpos)
}
func (br *bamRecord) isize() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.core.isize)
}
func (br *bamRecord) setIsize(isize int32) {
	if br.b == nil {
		panic(valueIsNil)
	}
	br.b.core.isize = C.int32_t(isize)
}
func (br *bamRecord) lAux() int32 {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int32(br.b.l_aux)
}
func (br *bamRecord) setLAux(lAux int32) {
	if br.b == nil {
		panic(valueIsNil)
	}
	br.b.l_aux = C.int(lAux)
}
func (br *bamRecord) dataLen() int {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int(br.b.data_len)
}
func (br *bamRecord) dataCap() int {
	if br.b == nil {
		panic(valueIsNil)
	}
	return int(br.b.m_data)
}
func (br *bamRecord) dataPtr() uintptr {
	if br.b == nil {
		panic(valueIsNil)
	}
	return uintptr(unsafe.Pointer(br.b.data))
}
func (br *bamRecord) dataUnsafe() []byte {
	if br.b == nil {
		panic(valueIsNil)
	}

	l := int(br.b.data_len)
	var data []byte
	sliceHeader := (*reflect.SliceHeader)(unsafe.Pointer(&data))
	sliceHeader.Cap = l
	sliceHeader.Len = l
	sliceHeader.Data = uintptr(unsafe.Pointer(br.b.data))

	return data
}
func (br *bamRecord) setDataUnsafe(data []byte) {
	if br.b == nil {
		panic(valueIsNil)
	}

	l := len(data)
	if br.dataCap() < len(data) {
		if br.b.data == nil {
			br.b.data = (*C.uint8_t)(unsafe.Pointer(C.malloc((C.size_t)(l))))
		} else {

			br.b.data = (*C.uint8_t)(unsafe.Pointer(C.realloc(unsafe.Pointer(br.b.data), (C.size_t)(l))))
		}
		if br.b.data == nil {
			panic(couldNotAllocate)
		}
	}

	var newData []byte
	sliceHeader := (*reflect.SliceHeader)(unsafe.Pointer(&newData))
	sliceHeader.Cap = l
	sliceHeader.Len = l
	sliceHeader.Data = uintptr(unsafe.Pointer(br.b.data))
	copy(newData, data)
}

// bamRecordFree C.free()s the contained bam1_t and its data, first checking for nil pointers.
func (br *bamRecord) bamRecordFree() {
	if br.b != nil {
		if br.b.data != nil {
			C.free(unsafe.Pointer(br.b.data))
		}
		C.free(unsafe.Pointer(br.b))
		br.b = nil
	}
}

// A samFile wraps a samfile_t.
type samFile struct {
	fp *C.samfile_t
}

// samOpen/samFdOpen open a SAM or BAM file with the given filename/fd, mode and optional auxilliary header.
// According to sam.h:
//
// mode matches /[rw](b?)(u?)(h?)([xX]?)/
//
//   'r' for reading,
//   'w' for writing,
//   'b' for BAM I/O,
//   'u' for uncompressed BAM output,
//   'h' for outputing header in SAM,
//   'x' for HEX flag and
//   'X' for string flag.
//
// If 'b' present, it must immediately follow 'r' or 'w'.
// Valid modes are "r", "w", "wh", "wx", "whx", "wX", "whX", "rb", "wb" and "wbu" exclusively.
//
// If mode[0] == 'w', aux must be a bamHeader.
// If mode[0] == 'r' && mode != "rb" and @SQ header lines in SAM are absent,
// aux must contain the name of a file listing of the reference sequences in SAM format.
// If neither of these conditions is met aux is not used.
func samOpen(filename, mode string, aux header) (sf *samFile, err error) {
	fn, m := C.CString(filename), C.CString(mode)
	defer C.free(unsafe.Pointer(fn))
	defer C.free(unsafe.Pointer(m))

	var auxAddr unsafe.Pointer
	switch a := aux.(type) {
	case textHeader:
		if len(a) > 0 {
			auxAddr = unsafe.Pointer(&a[0])
		} else {
			auxAddr = nil
		}
	case stringHeader:
		auxAddr = unsafe.Pointer(C.CString(string(a)))
		defer C.free(unsafe.Pointer(auxAddr))
	case *bamHeader:
		auxAddr = unsafe.Pointer(a.bh)
	default:
		if aux == nil {
			break
		}
		panic(fmt.Sprintf("boom: wrong type %T", aux))
	}

	fp, err := C.samopen(
		(*C.char)(unsafe.Pointer(fn)),
		(*C.char)(unsafe.Pointer(m)),
		unsafe.Pointer(auxAddr),
	)
	sf = &samFile{fp: (*C.samfile_t)(unsafe.Pointer(fp))}
	runtime.SetFinalizer(sf, (*samFile).samClose)

	return
}
func samFdOpen(fd uintptr, mode string, aux header) (sf *samFile, err error) {
	m := C.CString(mode)
	defer C.free(unsafe.Pointer(m))

	var auxAddr unsafe.Pointer
	switch a := aux.(type) {
	case textHeader:
		if len(a) > 0 {
			auxAddr = unsafe.Pointer(&a[0])
		} else {
			auxAddr = nil
		}
	case stringHeader:
		auxAddr = unsafe.Pointer(C.CString(string(a)))
		defer C.free(unsafe.Pointer(auxAddr))
	case *bamHeader:
		auxAddr = unsafe.Pointer(a.bh)
	default:
		if aux == nil {
			break
		}
		panic(fmt.Sprintf("boom: wrong type %T", aux))
	}

	fp, err := C.samdopen(
		C.int(fd),
		(*C.char)(unsafe.Pointer(m)),
		auxAddr,
	)
	sf = &samFile{fp: (*C.samfile_t)(unsafe.Pointer(fp))}
	runtime.SetFinalizer(sf, (*samFile).samClose)

	return
}

type bamTypeFlags int

const (
	// TODO: Curent definitions are a bit haphazard due to the underlying libbam defs. When tests exist, use 1<<iota.
	bamFile  bamTypeFlags = iota + 1         // File is a BAM file. Defined in sam.c TYPE_BAM
	readFile                                 // File is opened for reading. Defined in sam.c TYPE_READ
	hexFlags bamTypeFlags = C.BAM_OFHEX << 2 // Flags are in string format.
	strFlags bamTypeFlags = C.BAM_OFSTR << 2 // Flags are in hex format.
)

// fileType returns the type of file wrapped by the samFile struct.
func (sf *samFile) fileType() bamTypeFlags {
	if sf.fp != nil {
		return bamTypeFlags(sf.fp._type)
	}
	panic(valueIsNil)
}

// header returns the bamHeader wrapping the bam_header_t associated with sf.fp
func (sf *samFile) header() *bamHeader {
	if sf.fp == nil {
		return nil
	}
	return &bamHeader{bh: sf.fp.header}
}

// samClose closes the samFile, freeing the C data allocations as part of C.samclose.
func (sf *samFile) samClose() error {
	if sf.fp == nil {
		return valueIsNil
	}
	runtime.SetFinalizer(sf, nil)

	if h := sf.header(); h.bh.hash != nil {
		C.bam_destroy_header_hash(
			(*C.bam_header_t)(unsafe.Pointer(h.bh)),
		)
	}

	C.samclose((*C.samfile_t)(unsafe.Pointer(sf.fp)))
	sf.fp = nil

	return nil
}

// samRead reads and returns the next BAM record returning the number of bytes read,
// a *bamRecord containing the record data and any error that occurred.
func (sf *samFile) samRead() (n int, br *bamRecord, err error) {
	if sf.fp == nil {
		return 0, nil, valueIsNil
	}

	br, err = newBamRecord(nil)
	if err != nil {
		return
	}

	cn, err := C.samread(
		(*C.samfile_t)(unsafe.Pointer(sf.fp)),
		(*C.bam1_t)(unsafe.Pointer(br.b)),
	)
	n = int(cn)
	if n < 0 {
		err = io.EOF
	}

	return
}

// samWrite writes a BAM record represented by br, returning the number of bytes written
// and any error that occurred.
func (sf *samFile) samWrite(br *bamRecord) (n int, err error) {
	if sf.fp == nil || br.b == nil {
		return 0, valueIsNil
	}

	return int(C.samwrite(
		(*C.samfile_t)(unsafe.Pointer(sf.fp)),
		(*C.bam1_t)(unsafe.Pointer(br.b)),
	)), nil
}

// A bamIndex wraps a bam_index_t.
type bamIndex struct {
	idx *C.bam_index_t
}

// bamIndexBuild builds a BAM index file, filename.bai, from a bam file, filename. It returns an
// integer value (currently defined as always 0) and any error that occured.
func bamIndexBuild(filename string) (ret int, err error) {
	fn := C.CString(filename)
	defer C.free(unsafe.Pointer(fn))

	r, err := C.bam_index_build(
		(*C.char)(unsafe.Pointer(fn)),
	)

	return int(r), err
}

// bamIndexLoad loads a BAM index, returning a *bamIndex and any error that occurred.
// The error should be checked as a non-nil bamIndex is returned independent of error conditions.
// The bamIndex is created setting a finaliser that C.free()s the contained bam_index_t.
func bamIndexLoad(filename string) (bi *bamIndex, err error) {
	fn := C.CString(filename)
	defer C.free(unsafe.Pointer(fn))

	ip, err := C.bam_index_load(
		(*C.char)(unsafe.Pointer(fn)),
	)
	bi = &bamIndex{idx: (*C.bam_index_t)(unsafe.Pointer(ip))}
	runtime.SetFinalizer(bi, (*bamIndex).bamIndexDestroy)

	return
}

// bamIndexDestroy C.free()s the contained bam_index_t and its data, first checking for nil pointers.
func (bi *bamIndex) bamIndexDestroy() (err error) {
	if bi.idx == nil {
		return valueIsNil
	}

	C.bam_index_destroy(
		(*C.bam_index_t)(unsafe.Pointer(bi.idx)),
	)

	return
}

// A bamFetchFn is called on each bamRecord found by bamFetch. The return value is used to indicate
// the iteration is complete.
type bamFetchFn func(*bamRecord) bool

// bamFetch calls fn on all BAM records within the interval [beg, end) of the reference sequence
// identified by tid. Note that beg >= 0 || beg = 0.
func (sf *samFile) bamFetch(bi *bamIndex, tid, beg, end int, fn bamFetchFn) (ret int, err error) {
	if sf.fp == nil || bi.idx == nil {
		return 0, valueIsNil
	}

	if sf.fileType()&bamFile == 0 {
		return 0, notBamFile
	}

	fp := *(*C.bamFile)(unsafe.Pointer(&sf.fp.x))
	iter := C.bam_iter_query(bi.idx, C.int(tid), C.int(beg), C.int(end))
	var br *bamRecord
	for {
		br, err = newBamRecord(nil)
		if err != nil {
			return
		}
		ret = int(C.bam_iter_read(fp, iter, br.b))
		if ret < 0 {
			break
		}
		if fn(br) {
			break
		}
	}
	C.bam_iter_destroy(iter)

	return
}

// A bamFetchCFn is called on each bam1_t found by bamFetchC and the unsafe.Pointer is passed as a
// pointer to a store of user data. The integer return value is ignored internally by bam_fetch,
// but is specified in the libbam headers.
type bamFetchCFn func(*C.bam1_t, unsafe.Pointer) C.int

// bamFetchC calls fn on all BAM records within the interval [beg, end) of the reference sequence
// identified by tid. Note that beg >= 0 || beg = 0. data is passed to fn.
func (sf *samFile) bamFetchC(bi *bamIndex, tid, beg, end int, data unsafe.Pointer, fn bamFetchCFn) (ret int, err error) {
	if sf.fp == nil || bi.idx == nil {
		return 0, valueIsNil
	}

	if sf.fileType()&bamFile == 0 {
		return 0, notBamFile
	}

	r := C.bam_fetch(
		*(*C.bamFile)(unsafe.Pointer(&sf.fp.x)),
		bi.idx,
		C.int(tid),
		C.int(beg),
		C.int(end),
		data,
		(*[0]byte)(unsafe.Pointer(&fn)),
	)

	return int(r), nil
}

// Type header defines types that can be passed to samOpen as a SAM header or header filename.
type header interface {
	header() // No-op for interface definition.
}

// A bamHeader wraps a bam_header_t.
type bamHeader struct {
	bh *C.bam_header_t
}

// bamGetTid return the target id for for a reference sequence target matching the string, name.
func (bh *bamHeader) bamGetTid(name string) int {
	if bh.bh == nil {
		panic(valueIsNil)
	}

	sn := C.CString(name)
	defer C.free(unsafe.Pointer(sn))

	C.bam_init_header_hash( // This is idempotent - checks against NULL in bam_aux.c
		(*C.bam_header_t)(unsafe.Pointer(bh.bh)),
	)
	tid := C.bam_get_tid(
		(*C.bam_header_t)(unsafe.Pointer(bh.bh)),
		(*C.char)(unsafe.Pointer(sn)),
	)

	return int(tid)
}

// nTargets returns the number of reference sequence targets described in the BAM header.
func (bh *bamHeader) nTargets() int32 {
	if bh.bh != nil {
		return int32(bh.bh.n_targets)
	}
	panic(valueIsNil)
}

// targetNames returns a slice of strings containing the names of the reference sequence
// targets described in the BAM header.
func (bh *bamHeader) targetNames() (n []string) {
	if bh.bh != nil {
		n = make([]string, bh.bh.n_targets)
		l := int(bh.bh.n_targets)
		var nPtrs []*C.char
		sh := (*reflect.SliceHeader)(unsafe.Pointer(&nPtrs))
		sh.Cap = l
		sh.Len = l
		sh.Data = uintptr(unsafe.Pointer(bh.bh.target_name))

		for i, p := range nPtrs {
			n[i] = C.GoString(p)
		}

		return
	}
	panic(valueIsNil)
}

// targetLengths returns a slice of uint32 containing the lengths of the reference sequence
// targets described in the BAM header.
func (bh *bamHeader) targetLengths() []uint32 {
	if bh.bh != nil {
		l := int(bh.bh.n_targets)
		var unsafeLengths []uint32
		sh := (*reflect.SliceHeader)(unsafe.Pointer(&unsafeLengths))
		sh.Cap = l
		sh.Len = l
		sh.Data = uintptr(unsafe.Pointer(bh.bh.target_len))

		return append([]uint32(nil), unsafeLengths...)
	}
	panic(valueIsNil)
}

// text returns a string containing the full unparsed BAM header.
func (bh *bamHeader) text() (t string) {
	if bh.bh != nil {
		return C.GoStringN(bh.bh.text, C.int(bh.bh.l_text))
	}
	panic(valueIsNil)
}

// header is a no-op function required to allow *bamHeader to satisfy the header interface.
func (bh *bamHeader) header() {}

// stringHeader is a string representation of a filename of a SAM header file.
type stringHeader string

// header is a no-op function required to allow stringHeader to satisfy the header interface.
func (sh stringHeader) header() {}

// textHeader is a []byte representation of a filename of a SAM header file.
type textHeader []byte

// header is a no-op function required to allow textHeader to satisfy the header interface.
func (th textHeader) header() {}

const (
	Paired       Flags = paired       // The read is paired in sequencing, no matter whether it is mapped in a pair.
	ProperPair   Flags = properPair   // The read is mapped in a proper pair.
	Unmapped     Flags = unmapped     // The read itself is unmapped; conflictive with BAM_FPROPER_PAIR.
	MateUnmapped Flags = mateUnmapped // The mate is unmapped.
	Reverse      Flags = reverse      // The read is mapped to the reverse strand.
	MateReverse  Flags = mateReverse  // The mate is mapped to the reverse strand.
	Read1        Flags = read1        // This is read1.
	Read2        Flags = read2        // This is read2.
	Secondary    Flags = secondary    // Not primary alignment.
	QCFail       Flags = qCFail       // QC failure.
	Duplicate    Flags = duplicate    // Optical or PCR duplicate.
)

const (
	paired       = C.BAM_FPAIRED
	properPair   = C.BAM_FPROPER_PAIR
	unmapped     = C.BAM_FUNMAP
	mateUnmapped = C.BAM_FMUNMAP
	reverse      = C.BAM_FREVERSE
	mateReverse  = C.BAM_FMREVERSE
	read1        = C.BAM_FREAD1
	read2        = C.BAM_FREAD2
	secondary    = C.BAM_FSECONDARY
	qCFail       = C.BAM_FQCFAIL
	duplicate    = C.BAM_FDUP
)

// A Flags represents a BAM record's alignment FLAG field.
type Flags uint32

// String representation of BAM alignment flags:
//  0x001 - p - Paired
//  0x002 - P - ProperPair
//  0x004 - u - Unmapped
//  0x008 - U - MateUnmapped
//  0x010 - r - Reverse
//  0x020 - R - MateReverse
//  0x040 - 1 - Read1
//  0x080 - 2 - Read2
//  0x100 - s - Secondary
//  0x200 - f - QCFail
//  0x400 - d - Duplicate
//
// Note that flag bits are represented high order to the right.
func (f Flags) String() string {
	// If 0x01 is unset, no assumptions can be made about 0x02, 0x08, 0x20, 0x40 and 0x80
	const pairedMask = ProperPair | MateUnmapped | MateReverse | MateReverse | Read1 | Read2
	if f&1 == 0 {
		f &^= pairedMask
	}

	const flags = "pPuUrR12sfd"

	b := make([]byte, len(flags))
	for i, c := range flags {
		if f&(1<<uint(i)) != 0 {
			b[i] = byte(c)
		} else {
			b[i] = '-'
		}
	}

	return string(b)
}
