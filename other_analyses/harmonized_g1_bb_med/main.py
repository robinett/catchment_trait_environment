from harmonize import harmonize

def main():
    # what is bb g1 we want to convert from?
    bb_g1 = 8
    # what is the temperature we are assuming (in K)?
    T = 298
    # what is the rel humidity we are assuming (fractional)?
    rel_hum = 0.8

    h = harmonize()
    med_g1 = h.get_med_g1_from_bb_g1(bb_g1,T,rel_hum)
    out_str = (
        'for BB g1 of {}, temp of {}K, and rel humidty of {} ' +
        'the harmonized medlyn g1 is {}'
    )

    print(
        out_str.format(
            bb_g1,T,rel_hum,med_g1
        )
    )

if __name__ == '__main__':
    main()
