package vsop87_test

import (
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/soniakeys/vsop87"
)

var body = map[string]int{
	"MERCURY":    vsop87.Mercury,
	"VENUS":      vsop87.Venus,
	"MARS":       vsop87.Mars,
	"JUPITER":    vsop87.Jupiter,
	"SATURN":     vsop87.Saturn,
	"URANUS":     vsop87.Uranus,
	"NEPTUNE":    vsop87.Neptune,
	"EARTH-MOON": vsop87.EarthMoon,
}

func TestChk(t *testing.T) {
	vpath := os.Getenv("VSOP87")
	if vpath == "" {
		t.Log("Environment variable VSOP87 not set.")
		t.Skip("VSOP87 must be path to VSOP87 files.")
	}
	d, err := ioutil.ReadFile(filepath.Join(vpath, "vsop87.chk"))
	if err != nil {
		t.Fatal(err)
	}
	em, err := vsop87.NewEllipticModel(vpath, 0, 2451545)
	if err != nil {
		t.Fatal(err)
	}
	var e vsop87.Elliptic
	lines := strings.Split(string(d), "\n")
	t.Log(len(lines), "lines")
	for ln := 0; ln+3 < len(lines); {
		line := lines[ln]
		ln++ // report 1 based line numbers in log
		if line[0:8] != " VSOP87 " {
			t.Logf("%q", line[0:8])
			break
		}
		bs := strings.TrimSpace(line[8:21])
		b, ok := body[bs]
		if !ok {
			t.Fatal("ln: %d: %s unrecognized body", ln, bs)
		}
		jd, err := strconv.ParseFloat(strings.TrimSpace(line[24:34]), 64)
		if err != nil {
			t.Fatalf("ln %d: %v", ln, err)
		}
		t.Log("VSOP87", bs, jd)
		em.Pos(jd, b, &e)
		f := strings.Fields(lines[ln])
		ln++
		if len(f) != 9 {
			t.Errorf("ln %d: expected a k q", ln)
			t.Fatal("         got",lines[ln])
		}
		if s := chkStr(e.A); s != f[1] {
			t.Fatalf("ln %d: expected a = %s, got %s", ln, s, f[1])
		}
		if s := chkStr(e.K); s != f[4] {
			t.Fatalf("ln %d: expected k = %s, got %s", ln, s, f[4])
		}
		if s := chkStr(e.Q); s != f[7] {
			t.Fatalf("ln %d: expected q = %s, got %s", ln, s, f[7])
		}
		f = strings.Fields(lines[ln])
		ln++
		if len(f) != 9 {
			t.Fatalf("ln %d: expected l h p got\n%s", ln, lines[ln])
		}
		if s := chkStr(e.L); s != f[1] {
			t.Fatalf("ln %d: expected l = %s, got %s", ln, s, f[1])
		}
		if s := chkStr(e.H); s != f[4] {
			t.Fatalf("ln %d: expected h = %s, got %s", ln, s, f[4])
		}
		if s := chkStr(e.P); s != f[7] {
			t.Fatalf("ln %d: expected p = %s, got %s", ln, s, f[7])
		}
		ln++ // skip blank line
	}
}

func chkStr(f float64) string {
	switch s := fmt.Sprintf("%.10f", f); {
	case s[0] == '0':
		return s[1:]
	case s[:2] == "-0":
		return "-" + s[2:]
	default:
		return s
	}
}
