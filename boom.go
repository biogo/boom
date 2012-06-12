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
*/
import "C"

import (
	"encoding/binary"
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

type bamRecord struct {
	b *C.bam1_t
}

func newBamRecord(b *C.bam1_t) (br *bamRecord) {
	br = &bamRecord{b}
	runtime.SetFinalizer(br, (*bamRecord).bamRecordFree)

	return
}

func (br *bamRecord) tid() int32 {
	if br.b != nil {
		return int32(br.b.core.tid)
	}
	panic(valueIsNil)
}
func (br *bamRecord) pos() int32 {
	if br.b != nil {
		return int32(br.b.core.pos)
	}
	panic(valueIsNil)
}
func (br *bamRecord) bin() uint {
	if br.b != nil {
		return uint(br.b.core.bin)
	}
	panic(valueIsNil)
}
func (br *bamRecord) qual() byte {
	if br.b != nil {
		return byte(br.b.core.qual)
	}
	panic(valueIsNil)
}
func (br *bamRecord) lQname() byte {
	if br.b != nil {
		return byte(br.b.core.l_qname)
	}
	panic(valueIsNil)
}
func (br *bamRecord) flag() Flags {
	if br.b != nil {
		return Flags(br.b.core.flag)
	}
	panic(valueIsNil)
}
func (br *bamRecord) nCigar() uint {
	if br.b != nil {
		return uint(br.b.core.n_cigar)
	}
	panic(valueIsNil)
}
func (br *bamRecord) lQseq() int32 {
	if br.b != nil {
		return int32(br.b.core.l_qseq)
	}
	panic(valueIsNil)
}
func (br *bamRecord) mtid() int32 {
	if br.b != nil {
		return int32(br.b.core.mtid)
	}
	panic(valueIsNil)
}
func (br *bamRecord) mpos() int32 {
	if br.b != nil {
		return int32(br.b.core.mpos)
	}
	panic(valueIsNil)
}
func (br *bamRecord) isize() int32 {
	if br.b != nil {
		return int32(br.b.core.isize)
	}
	panic(valueIsNil)
}
func (br *bamRecord) lAux() int {
	if br.b != nil {
		return int(br.b.l_aux)
	}
	panic(valueIsNil)
}
func (br *bamRecord) dataLen() int {
	if br.b != nil {
		return int(br.b.data_len)
	}
	panic(valueIsNil)
}
func (br *bamRecord) dataCap() int {
	if br.b != nil {
		return int(br.b.m_data)
	}
	panic(valueIsNil)
}
func (br *bamRecord) dataPtr() uintptr {
	if br.b != nil {
		return uintptr(unsafe.Pointer(br.b.data))
	}
	panic(valueIsNil)
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
func (br *bamRecord) bamRecordFree() {
	if br.b != nil {
		C.free(unsafe.Pointer(br.b))
		br.b = nil
	}
}

type samFile struct {
	fp *C.samfile_t
}

func (sf *samFile) fileType() int {
	if sf.fp != nil {
		return int(sf.fp._type)
	}
	panic(valueIsNil)
}

func samOpen(filename, mode string, aux header) (sf *samFile, err error) {
	fn, m := C.CString(filename), C.CString(mode)
	defer C.free(unsafe.Pointer(fn))
	defer C.free(unsafe.Pointer(m))

	var auxAddr uintptr
	switch aux.(type) {
	case textHeader:
		xv := reflect.ValueOf(aux)
		if xv.Len() > 0 {
			auxAddr = xv.Index(0).UnsafeAddr()
		} else {
			auxAddr = 0
		}
	case *bamHeader:
		auxAddr = reflect.ValueOf(aux).UnsafeAddr()
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

func (sf *samFile) header() *bamHeader {
	return &bamHeader{bh: sf.fp.header}
}

func (sf *samFile) samClose() error {
	if sf.fp == nil {
		return valueIsNil
	}

	C.samclose((*C.samfile_t)(unsafe.Pointer(sf.fp)))

	return nil
}

func (sf *samFile) samRead() (n int, br *bamRecord, err error) {
	if sf.fp == nil {
		return 0, nil, valueIsNil
	}

	b := (*C.bam1_t)(unsafe.Pointer(C.malloc((C.size_t)(unsafe.Sizeof(C.bam1_t{})))))

	if b == nil {
		return 0, nil, couldNotAllocate
	}
	*b = C.bam1_t{}
	br = newBamRecord(b)

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

func (sf *samFile) samWrite(br *bamRecord) (n int, err error) {
	if sf.fp == nil || br.b == nil {
		return 0, valueIsNil
	}

	return int(C.samwrite(
		(*C.samfile_t)(unsafe.Pointer(sf.fp)),
		(*C.bam1_t)(unsafe.Pointer(br.b)),
	)), nil
}

type bamIndex struct {
	idx *C.bam_index_t
}

func bamIndexBuild(filename string) (ret int, err error) {
	fn := C.CString(filename)
	defer C.free(unsafe.Pointer(fn))

	r, err := C.bam_index_build(
		(*C.char)(unsafe.Pointer(fn)),
	)

	return int(r), err
}

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

func (bi *bamIndex) bamIndexDestroy() (err error) {
	if bi.idx == nil {
		return valueIsNil
	}

	C.bam_index_destroy(
		(*C.bam_index_t)(unsafe.Pointer(bi.idx)),
	)

	return
}

type BamFetchFn func(*bamRecord, interface{}) int

func (sf *samFile) bamFetch(bi *bamIndex, tid, beg, end int, data interface{}, fn BamFetchFn) (ret int, err error) {
	if sf.fp == nil || bi.idx == nil {
		return 0, valueIsNil
	}

	if sf.fileType()&1 != 1 {
		return 0, notBamFile
	}

	dv := reflect.ValueOf(data)
	if dv.CanAddr() {
		return 0, cannotAddr
	}

	f := func(b *C.bam1_t, _ *[0]byte) C.int {
		return C.int(fn(&bamRecord{(*C.bam1_t)(unsafe.Pointer(b))}, data))
	}

	r := C.bam_fetch(
		(*C.BGZF)(unsafe.Pointer(sf.fp)),
		bi.idx,
		C.int(tid),
		C.int(beg),
		C.int(end),
		unsafe.Pointer(dv.UnsafeAddr()),
		(*[0]byte)(unsafe.Pointer(&f)),
	)

	return int(r), nil
}

type header interface {
	header()
}

type bamHeader struct {
	bh *C.bam_header_t
}

func (bh *bamHeader) bamGetTid(name string) (int, error) {
	sn := C.CString(name)
	defer C.free(unsafe.Pointer(sn))

	tid, err := C.bam_get_tid(
		(*C.bam_header_t)(unsafe.Pointer(bh.bh)),
		(*C.char)(unsafe.Pointer(sn)),
	)

	return int(tid), err
}

func (bh *bamHeader) nTargets() int32 {
	if bh.bh != nil {
		return int32(bh.bh.n_targets)
	}
	panic(valueIsNil)
}

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

func (bh *bamHeader) text() (t string) {
	if bh.bh != nil {
		return C.GoStringN(bh.bh.text, C.int(bh.bh.l_text))
	}
	panic(valueIsNil)
}

func (bh *bamHeader) header() {}

type textHeader []byte

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
	paired             = C.BAM_FPAIRED
	properPair         = C.BAM_FPROPER_PAIR
	unmapped           = C.BAM_FUNMAP
	mateUnmapped       = C.BAM_FMUNMAP
	reverse            = C.BAM_FREVERSE
	mateReverse        = C.BAM_FMREVERSE
	read1              = C.BAM_FREAD1
	read2              = C.BAM_FREAD2
	secondary          = C.BAM_FSECONDARY
	qCFail             = C.BAM_FQCFAIL
	duplicate          = C.BAM_FDUP
)

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
