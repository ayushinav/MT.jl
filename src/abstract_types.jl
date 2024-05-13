
"""
Abstract model type that is the supertype of all `model`s in the package.
"""
abstract type AbstractModel end

"""
Abstract model type that is the supertype of all geophysical models in the package.
"""
abstract type AbstractGeophyModel <: AbstractModel end

"""
Abstract model type that is the supertype of all `model`s in the package.
"""
abstract type AbstractResponse end

"""
Abstract model type that is the supertype of all geophysical models in the package.
"""
abstract type AbstractGeophyResponse <: AbstractResponse end