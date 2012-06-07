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
*/
import "C"

import (
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
func (br *bamRecord) flag() uint {
	if br.b != nil {
		return uint(br.b.core.flag)
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
func (br *bamRecord) data() uintptr {
	if br.b != nil {
		return uintptr(unsafe.Pointer(br.b.data))
	}
	panic(valueIsNil)
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

func (bh *bamHeader) header() {}

type textHeader []byte

func (th textHeader) header() {}
