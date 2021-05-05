def square_residue(elem, ideal):
    if elem.parent() != ZZ:
        elem = ideal.reduce(elem)
        if elem.is_zero():
            return True

        for res in ideal.residues():
            if ideal.reduce(res^2) == elem:
                return True

        return False
    else:
        elem = elem % ideal
        if elem == 0:
            return True

        for i in range(ideal):
            if i^2 % ideal == elem:
                return True

        return False


def squarefree_residue(elem, ideal):
    ideal_factors = ideal.factor()
    ideal_contains_square_factors = False
    for fctr in ideal_factors:
        if fctr[1] > 1:
            ideal_contains_square_factors = True
            break

    if not ideal_contains_square_factors:
        return True

    if elem.is_zero():
        return False

    if elem.parent() == ZZ:
        elem_factors = elem.factor()
    else:
        elem_factors = elem.parent().ideal(elem).factor()
    for elem_fctr in elem_factors:
        for fctr in ideal_factors:
            if fctr[1] == 1:
                continue
            if elem_fctr[0] == fctr[0] and elem_fctr[1] > 1:
                return False

    return True


def factor_mult(ideal):
    sum = 0
    for fctr in ideal.factor():
        sum += fctr[1]
    return sum


def ideals_dividing(ideal):
    factors = ideal.factor()
    if type(ideal) == Integer:
        ideal_list = [1]
    else:
        ideal_list = [ideal.number_field().ideal(1)]
    for fctr in factors:
        temp_lis = deepcopy(ideal_list)
        for i in range(fctr[1]):
            for idl in temp_lis:
                ideal_list.append(idl * (fctr[0]^(i+1)))

    ideal_list.sort(key=factor_mult, reverse=True)
    return ideal_list


def two_part_disc(elem):
    parent_ring = elem.parent()
    if parent_ring == ZZ:
        two = 2
        if elem%4 == 0:
            elem = elem//4
        elem_ideal = elem
    else:
        elem_ideal = parent_ring.ideal(elem)
        two = parent_ring.ideal(2)

    if not squarefree_residue(elem, two^2):
        return None

    two_part = two^2

    two_divs = ideals_dividing(two)

    for fctr in elem_ideal.factor():
        if fctr[0].divides(two):
            two_part *= fctr[0]
            two_divs = list(filter(lambda x: not fctr[0].divides(x), two_divs))

    for idl in two_divs:
        if square_residue(elem, idl^2):
            return two_part/(idl^2)


def local_quad_reps(K, ref_elem = None):
    rep_list = []
    cumulative_ideal = K.ideal(1)
    for p, e in K.ideal(2).factor():
        pi = K.uniformizer(p)
        p_rep_list = []
        p_squares = []
        temp_id = p^(2*e+1)
        p_elem_list = [temp_id.reduce(elem) for elem in temp_id.invertible_residues()]
        for i in range(len(p_elem_list)):
            p_elem_list[i] = temp_id.reduce(p_elem_list[i])
            if temp_id.reduce(p_elem_list[i]^2) not in p_squares:
                p_squares.append(temp_id.reduce(p_elem_list[i]^2))

        if ref_elem is None:
            for p_elem in p_elem_list:
                for square in p_squares:
                    if temp_id.reduce(p_elem * square) in p_rep_list:
                        break
                else:
                    p_rep_list.append(p_elem)
                    p_rep_list.append((temp_id*p).reduce(pi * p_elem))
        else:
            try:
                for p_elem in p_elem_list:
                    for square in p_squares:
                        if temp_id.reduce(p_elem * square) == temp_id.reduce(ref_elem):
                            p_rep_list.append(p_elem)
                            raise Exception
                        elif (temp_id*p).reduce(pi*p_elem*square) == \
                            (temp_id*p).reduce(ref_elem):
                            p_rep_list.append(pi*p_elem)
                            raise Exception
            except:
                pass

        if len(rep_list) == 0:
            rep_list = p_rep_list
        else:
            temp_list = []
            for r_elem in rep_list:
                for p_elem in p_rep_list:
                    temp_list.append(
                        (cumulative_ideal * temp_id*p).reduce(
                            K.solve_CRT([r_elem, p_elem], [cumulative_ideal, temp_id*p])
                         )
                    )
            rep_list = temp_list
        cumulative_ideal *= temp_id * p

    return rep_list


###CODE BELOW HERE GENERATES TABLE OF LOCAL WEIGHTS
for i in [2, -2, 5, -5, 10, -10, -1, -7]:
    K = QuadraticField(i)
    disc_weight = 0
    D_L_K_rows = [
        [0,0,0],
        [0,0,0],
        [0,0,0],
        [0,0,0],
        [0,0,0],
        [0,0,0]
    ]
    for rep in local_quad_reps(K):
        relative_disc_2 = two_part_disc(rep)
        disc_weight += 1/relative_disc_2.norm()
        if relative_disc_2.norm() == 1:
            two_power = 0
        else:
            two_power = relative_disc_2.norm().factor()[0][1]
        flipped_disc_2 = two_part_disc(rep.norm()%16)
        if flipped_disc_2.norm() == 1:
            flipped_two_power = 0
        else:
            flipped_two_power = flipped_disc_2.norm().factor()[0][1]

        if two_power == 0:
            row_index = 0
        else:
            row_index = two_power - 1
        if flipped_two_power == 0:
            col_index = 0
        else:
            col_index = flipped_two_power - 1

        D_L_K_rows[row_index][col_index] += 1/relative_disc_2.norm()

    for x in range(6):
        for y in range(3):
             D_L_K_rows[x][y] /= disc_weight

    print(K)
    table(D_L_K_rows)
