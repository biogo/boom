// Copyright ©2012 The bíogo.boom Authors. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package boom

// BuildIndex builds a BAM index file, filename.bai, from a sorted BAM file, filename.
// It returns any error that occured.
func BuildIndex(file string) (err error) {
	_, err = bamIndexBuild(file)
	return
}

// An Index represents an in memory BAM index.
type Index struct {
	*bamIndex
}

// LoadIndex loads a BAM index file, and returns the index in i if no error occurred.
// If an error occurred i is returned nil and the error is returned.
func LoadIndex(file string) (i *Index, err error) {
	bi, err := bamIndexLoad(file)
	if err == nil {
		i = &Index{bi}
	}
	return
}
