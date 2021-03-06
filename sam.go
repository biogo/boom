// Copyright ©2012 The bíogo Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package boom

import (
	"os"
)

// A SAMFile represents a SAM (text Sequence Alignment/Map) file.
type SAMFile struct {
	*samFile
}

var tWModes = [2]string{"w", "wh"}

// OpenSAMFile opens the file, f as a SAM file.
// If an error occurrs it is returned with a nil SAMFile pointer.
// The valid values of mode and ref are described in the overview and are
// derived from the samtools documentation.
func OpenSAMFile(f *os.File, mode string, ref *Header) (b *SAMFile, err error) {
	var h header
	if ref != nil {
		h = ref.bamHeader
	}
	sf, err := samFdOpen(f.Fd(), mode, h)
	if err != nil {
		return
	}
	return &SAMFile{sf}, nil
}

// OpenSAM opens the file, filename as a SAM file.
// If an error occurrs it is returned with a nil SAMFile pointer.
func OpenSAM(filename, ref string) (s *SAMFile, err error) {
	h := textHeader(ref)
	sf, err := samOpen(filename, "r", h)
	if err != nil {
		return
	}
	return &SAMFile{sf}, nil
}

// CreateBAM opens a file, filename for writing. ref is required to point to a valid Header.
func CreateSAM(filename string, ref *Header, dh bool) (s *SAMFile, err error) {
	var mode string
	if dh {
		mode = tWModes[1]
	} else {
		mode = tWModes[0]
	}
	if ref == nil {
		return nil, noHeader
	}
	sf, err := samOpen(filename, mode, ref.bamHeader)
	if err != nil {
		return
	}
	return &SAMFile{sf}, nil
}

// Close closes the SAMFile, freeing any associated data.
func (self *SAMFile) Close() error {
	if self == nil {
		return nil
	}
	return self.samClose()
}

// Read reads a single SAM record and returns this or any error, and the number of bytes read.
func (self *SAMFile) Read() (r *Record, n int, err error) {
	n, br, err := self.samRead()
	r = &Record{bamRecord: br, marshalled: true}
	return
}

// Write writes a BAM record, r, returning the number of bytes written and any error that occurred.
func (self *SAMFile) Write(r *Record) (n int, err error) {
	if r.marshalled == false {
		r.setDataUnsafe(r.marshalData())
		r.marshalled = true
	}
	return self.samWrite(r.bamRecord)
}

// RefID returns the tid corresponding to the string chr and true if a match is present.
// If no matching tid is found -1 and false are returned.
func (self *SAMFile) RefID(chr string) (id int, ok bool) {
	id = self.header().bamGetTid(chr)
	if id < 0 {
		return
	}
	ok = true

	return
}

// Header returns a pointer to the BAM file's header.
func (self *SAMFile) Header() *Header {
	return &Header{self.header()}
}

// Targets returns the number of reference sequences described in the SAMFile's header.
func (self *SAMFile) Targets() int {
	h := self.header()
	if h == nil {
		return -1
	}
	return int(h.nTargets())
}

// RefNames returns a slice of strings containing the names of reference sequences described
// in the SAM file's header.
func (self *SAMFile) RefNames() []string {
	h := self.header()
	if h == nil {
		return nil
	}
	return h.targetNames()
}

// RefLengths returns a slice of integers containing the lengths of reference sequences described
// in the SAM file's header.
func (self *SAMFile) RefLengths() []uint32 {
	h := self.header()
	if h == nil {
		return nil
	}
	return h.targetLengths()
}

// Text returns the unparsed text of the SAM header as a string.
func (self *SAMFile) Text() string {
	h := self.header()
	if h == nil {
		return ""
	}
	return h.text()
}
