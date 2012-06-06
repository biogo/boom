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
		r, n, err := sf.Read()
		if err != nil {
			break
		}
		fmt.Println(n, r.Start(), r.End(), r.Len(), r.Quality(), r.Flags())
	}
}
