# Codes of the RTDR envelope piece types as .dump_rtdr_envelope_cpp returns
# them: the positions in the PieceType enum of src/mhn_rtdr.h, counted from 0.
PIECE_PLATEAU <- 0L            # constant at the density's peak value
PIECE_EXP_RIGHT <- 1L          # exponential right tail [a, Inf), anchored at a
PIECE_EXP_LEFT <- 2L           # exponential left tail (-Inf, a], anchored at a
PIECE_EXP_BOUNDED <- 3L        # exponential piece on [a, b], anchored at a
PIECE_SECANT <- 4L             # secant of log g on [a, b], anchored at a
PIECE_TNEGHALF_LEFT <- 5L      # inverse-square tangent at t = a on (-Inf, b - aux]
PIECE_TNEGHALF_RIGHT <- 6L     # inverse-square tangent at t = a on [b + aux, Inf)
PIECE_EXP_BOUNDED_RIGHT <- 7L  # exponential piece anchored at its right end
