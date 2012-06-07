package main

import (
	"fmt"
	"github.com/kortschak/boom"
)

func main() {
	sf, err := boom.OpenSAM("test.sam", "")
	if err != nil {
		panic(err)
	}
	for {
		r, _, err := sf.Read()
		if err != nil {
			break
		}
		fmt.Println(r.Start(), r.End(), r.Len(), r.Strand(), r.Quality(), r.Flags())
	}
}
