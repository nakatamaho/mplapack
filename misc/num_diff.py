import sys
import argparse
from decimal import Decimal, getcontext, InvalidOperation

# Set precision high enough for DD (31 digits) and QD (62 digits)
getcontext().prec = 70

def is_decimal(value):
    try:
        Decimal(value)
        return True
    except (InvalidOperation, ValueError):
        return False

def compare_files(file_ref, file_test, rel_tol, abs_tol):
    success = True
    rel_tol = Decimal(str(rel_tol))
    abs_tol = Decimal(str(abs_tol))
    
    with open(file_ref, 'r') as f_ref, open(file_test, 'r') as f_test:
        for line_no, (line1, line2) in enumerate(zip(f_ref, f_test), 1):
            tokens1 = line1.split()
            tokens2 = line2.split()

            if len(tokens1) != len(tokens2):
                print(f"Error at line {line_no}: Column count mismatch.")
                return False

            for col_no, (t1, t2) in enumerate(zip(tokens1, tokens2), 1):
                if t1 == t2:
                    continue
                
                if is_decimal(t1) and is_decimal(t2):
                    v1, v2 = Decimal(t1), Decimal(t2)
                    diff = abs(v1 - v2)
                    
                    # Calculate relative error
                    max_v = max(abs(v1), abs(v2))
                    rel_err = diff / max_v if max_v != 0 else diff
                    
                    if rel_err > rel_tol and diff > abs_tol:
                        print(f"Numerical Mismatch at Line {line_no}, Col {col_no}:")
                        print(f"  Expected: {t1}")
                        print(f"  Actual:   {t2}")
                        print(f"  RelErr:   {rel_err:.2e}")
                        success = False
                else:
                    print(f"String Mismatch at Line {line_no}, Col {col_no}:")
                    print(f"  Expected: '{t1}'")
                    print(f"  Actual:   '{t2}'")
                    success = False
            
            if not success:
                return False
                
    return True

def main():
    parser = argparse.ArgumentParser(description="High-precision numerical diff tool.")
    parser.add_argument("file_ref", help="Reference output file")
    parser.add_argument("file_test", help="Test output file")
    parser.add_argument("--tol", type=str, default="1e-30", help="Relative tolerance")
    parser.add_argument("--abstol", type=str, default="0", help="Absolute tolerance")

    args = parser.parse_args()

    # Pass tolerance as string to Decimal to maintain precision
    if compare_files(args.file_ref, args.file_test, args.tol, args.abstol):
        print("PASS: Files match within specified tolerance.")
        sys.exit(0)
    else:
        print("FAIL: Significant differences detected.")
        sys.exit(1)

if __name__ == "__main__":
    main()
