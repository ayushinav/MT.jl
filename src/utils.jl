# utils to help bump for Abstract types

import Base:zero
zero(x::response) = response([zero(getfield(x, k)) for k in fieldnames(response)]...)
