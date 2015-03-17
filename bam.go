// Copyright ©2012 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package boom

import (
	"os"
)

// A BAMFile represents a BAM (Binary Sequence Alignment/Map) file.
type BAMFile struct {
	*samFile
}

var bWModes = [2]string{"wb", "wbu"}

// OpenBAMFile opens the file, f as a BAM file.
// If an error occurrs it is returned with a nil BAMFile pointer.
// The valid values of mode and ref are described in the overview and are derived
// from the samtools documentation.
func OpenBAMFile(f *os.File, mode string, ref *Header) (b *BAMFile, err error) {
	var h header
	if ref != nil {
		h = ref.bamHeader
	}
	sf, err := samFdOpen(f.Fd(), mode, h)
	if err != nil {
		return
	}
	return &BAMFile{sf}, nil
}

// OpenBAM opens the file, filename as a BAM file.
// If an error occurrs it is returned with a nil BAMFile pointer.
func OpenBAM(filename string) (b *BAMFile, err error) {
	sf, err := samOpen(filename, "rb", nil)
	if err != nil {
		return
	}
	return &BAMFile{sf}, nil
}

// CreateBAM opens a file, filename for writing. ref is required to point to a valid Header.
// If comp is true, compression is used.
func CreateBAM(filename string, ref *Header, comp bool) (b *BAMFile, err error) {
	var mode string
	if comp {
		mode = bWModes[0]
	} else {
		mode = bWModes[1]
	}
	if ref == nil {
		return nil, noHeader
	}
	sf, err := samOpen(filename, mode, ref.bamHeader)
	if err != nil {
		return
	}
	return &BAMFile{sf}, nil
}

// Close closes the BAMFile, freeing any associated data.
func (self *BAMFile) Close() error {
	if self == nil {
		return nil
	}
	return self.samClose()
}

// Read reads a single BAM record and returns this or any error, and the number of bytes read.
func (self *BAMFile) Read() (r *Record, n int, err error) {
	n, br, err := self.samRead()
	r = &Record{bamRecord: br, marshalled: true}
	return
}

// Write writes a BAM record, r, returning the number of bytes written and any error that occurred.
func (self *BAMFile) Write(r *Record) (n int, err error) {
	if r.marshalled == false {
		r.setDataUnsafe(r.marshalData())
		r.marshalled = true
	}
	return self.samWrite(r.bamRecord)
}

// RefID returns the tid corresponding to the string chr and true if a match is present.
// If no matching tid is found -1 and false are returned.
func (self *BAMFile) RefID(chr string) (id int, ok bool) {
	id = self.header().bamGetTid(chr)
	if id < 0 {
		return
	}
	ok = true

	return
}

// Header returns a pointer to the BAM file's header.
func (self *BAMFile) Header() *Header {
	return &Header{self.header()}
}

// Targets returns the number of reference sequences described in the SAMFile's header.
func (self *BAMFile) Targets() int {
	h := self.header()
	if h == nil {
		return -1
	}
	return int(h.nTargets())
}

// RefNames returns a slice of strings containing the names of reference sequences described
// in the BAM file's header.
func (self *BAMFile) RefNames() []string {
	h := self.header()
	if h == nil {
		return nil
	}
	return h.targetNames()
}

// RefLengths returns a slice of integers containing the lengths of reference sequences described
// in the BAM file's header.
func (self *BAMFile) RefLengths() []uint32 {
	h := self.header()
	if h == nil {
		return nil
	}
	return h.targetLengths()
}

// Text returns the unparsed text of the BAM header as a string.
func (self *BAMFile) Text() string {
	h := self.header()
	if h == nil {
		return ""
	}
	return h.text()
}

// A FetchFn is called on each Record found by Fetch. Returning a true done value breaks from the
// iterator.
type FetchFn func(*Record) (done bool)

// Fetch calls fn on all BAM records within the interval [beg, end) of the reference sequence
// identified by chr. Note that beg >= 0 || beg = 0. The Record value passed by pointer to fn is reused
// each iteration and is unusable after Fetch returns, so the values should not be stored.
func (self *BAMFile) Fetch(i *Index, tid int, beg, end int, fn FetchFn) (ret int, err error) {
	f := func(b *bamRecord) bool {
		return fn(&Record{bamRecord: b})
	}

	return self.bamFetch(i.bamIndex, tid, beg, end, f)
}
