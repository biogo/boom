package main

import (
	"fmt"
	"github.com/kortschak/boom"
)

func main() {
	bf, err := boom.OpenBAM("test.bam")
	if err != nil {
		panic(err)
	}
	fmt.Printf("%s %d\n%s", sf.ReferenceNames(), sf.ReferenceLengths(), sf.Text())
	for {
		r, _, err := bf.Read()
		if err != nil {
			break
		}
		fmt.Printf("%d %s %d-%d (%d) %v %d %d %s %v Fl:%v %q\n", r.ReferenceID(), r.Name(), r.Start(), r.End(), r.Len(), r.Cigar(), r.Strand(), r.Score(), r.Seq(), r.Quality(), r.Flags(), r.Tags())
	}
}
