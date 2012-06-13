package main

import (
	"fmt"
	"code.google.com/p/biogo.boom"
)

func main() {
	boom.Verbosity(0)
	sf, err := boom.OpenSAM("test.sam", "")
	if err != nil {
		panic(err)
	}
	fmt.Printf("%s %d\n%s", sf.ReferenceNames(), sf.ReferenceLengths(), sf.Text())
	for {
		r, _, err := sf.Read()
		if err != nil {
			break
		}
		fmt.Printf("%d %s %d-%d (%d) %v %d %d %s %v Fl:%v %q\n", r.ReferenceID(), r.Name(), r.Start(), r.End(), r.Len(), r.Cigar(), r.Strand(), r.Score(), r.Seq(), r.Quality(), r.Flags(), r.Tags())
	}
}
